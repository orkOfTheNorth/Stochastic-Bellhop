"""
Stochastic Bellhop UQ — Interactive Streamlit viewer.

Run:  streamlit run ui/app.py  (from stochastic_TL_maps_project/)
      or
      streamlit run app.py     (from ui/)

Layout
------
Sidebar   : scenario, param, dist, FOM, detection threshold, overlay method
Left col  : EX(TL) heatmap with detection-probability overlay
Right col : Q-TIP point inspector — stats table + MC histogram
             Point selected via Range / Depth sliders below the map.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import scipy.io
import streamlit as st
from scipy.special import ndtr as norm_cdf  # Φ(x)

# ── Path resolution ────────────────────────────────────────────────────────────
# Works whether launched from ui/ or from stochastic_TL_maps_project/
HERE = Path(__file__).parent
PROJ_DIR = HERE if (HERE / "config.json").exists() else HERE.parent

# ── Streamlit page config ──────────────────────────────────────────────────────
st.set_page_config(
    page_title="Stochastic Bellhop UQ",
    layout="wide",
    initial_sidebar_state="expanded",
)


# ══════════════════════════════════════════════════════════════════════════════
# Config helpers
# ══════════════════════════════════════════════════════════════════════════════

@st.cache_resource
def load_config() -> dict:
    with open(PROJ_DIR / "config.json") as f:
        return json.load(f)


def get_fom(cfg: dict, scen: str) -> float:
    for s in cfg["scenarios"]:
        if s["name"] == scen:
            return float(s.get("FOM_dB", cfg["nominal"]["FOM_dB"]))
    return float(cfg["nominal"]["FOM_dB"])


PARAMS     = ["freq", "zS", "svp"]
PARAM_LBLS = {
    "freq": "Frequency (Hz)",
    "zS":   "Source depth (m)",
    "svp":  "SVP shift",
}


# ══════════════════════════════════════════════════════════════════════════════
# Data loaders (all cached)
# ══════════════════════════════════════════════════════════════════════════════

def _mat(path: Path, *keys):
    """Load a .mat file; return dict with only requested keys (or all if none given)."""
    if not path.exists():
        return None
    d = scipy.io.loadmat(str(path), squeeze_me=True)
    if keys:
        return {k: d[k] for k in keys if k in d}
    return d


@st.cache_data(show_spinner=False)
def load_mc(scen: str, dist: str, param: str) -> dict | None:
    p = PROJ_DIR / "Methods" / "MC" / scen / dist / "results" / f"MC_{param}.mat"
    d = _mat(p)
    if d is None:
        return None
    r = np.atleast_1d(d["r_km"]).ravel().astype(float)
    z = np.atleast_1d(d["z_m"]).ravel().astype(float)
    EX  = np.array(d["MC_EX"],  dtype=float)
    Var = np.array(d["MC_Var"], dtype=float)
    Pd  = d.get("MC_P_detect")
    P_detect = np.array(Pd, dtype=float) if Pd is not None else np.zeros_like(EX)
    return {"EX": EX, "Var": Var, "P_detect": P_detect, "r_km": r, "z_m": z}


@st.cache_data(show_spinner=False)
def load_delta(scen: str, dist: str, param: str) -> dict | None:
    """Load single-parameter delta results (delta_<param>.mat)."""
    p = PROJ_DIR / "Methods" / "Delta" / scen / dist / "results" / f"delta_{param}.mat"
    d = _mat(p)
    if d is None:
        return None
    return {
        "TL_expected":  np.array(d.get("TL_expected",  np.zeros((1,1))), dtype=float),
        "Vs_1st":       np.array(d.get("Vs_1st",       np.zeros((1,1))), dtype=float),
        "Var_TL":       np.array(d.get("Var_TL",       np.zeros((1,1))), dtype=float),
        "DeltaLN3_prob":np.array(d.get("DeltaLN3_prob",np.zeros((1,1))), dtype=float),
        "Cheb_lb_s":    np.array(d.get("Cheb_lb_s",    np.zeros((1,1))), dtype=float),
    }


@st.cache_data(show_spinner=False)
def load_pce(scen: str, dist: str) -> dict | None:
    """Load PCE results for all 3 params at once."""
    p = PROJ_DIR / "Methods" / "PCE" / scen / dist / "results" / "pce_results.mat"
    d = _mat(p)
    if d is None:
        return None
    Pd = d.get("Pdetect_best")
    if Pd is None:
        return None
    # MATLAB {1×3} cell → numpy object array shape (3,) after squeeze_me=True
    Pd = np.asarray(Pd)
    pdet = []
    for pi in range(3):
        try:
            arr = np.array(Pd.flat[pi], dtype=float)
        except Exception:
            arr = np.zeros((1, 1))
        pdet.append(arr)
    kstar_raw = d.get("Kstar")
    kstar = (np.atleast_1d(kstar_raw).ravel().astype(int)
             if kstar_raw is not None else np.zeros(3, int))
    return {"Pdetect": pdet, "Kstar": kstar}   # pdet[pi] → [Nz×Nr]


@st.cache_data(show_spinner=False)
def load_tl_cache(scen: str, param: str) -> np.ndarray | None:
    """Load raw 300-sample TL cube from base-distribution cache."""
    base = "Normal_10pct"
    candidates = [
        PROJ_DIR / "Cache" / scen / base / f"TL_{param}.mat",
        PROJ_DIR / "Cache" / scen / base / f"TL_{param}_N300.mat",
    ]
    for p in candidates:
        d = _mat(p)
        if d is not None:
            key = "TL_save" if "TL_save" in d else ("TL_all" if "TL_all" in d else None)
            if key:
                return np.array(d[key], dtype=float)   # [Nz × Nr × N]
    return None


# ══════════════════════════════════════════════════════════════════════════════
# Sidebar
# ══════════════════════════════════════════════════════════════════════════════

cfg = load_config()

SCENARIOS = [s["name"] for s in cfg["scenarios"]]
DISTS     = [d["name"] for d in cfg["distributions"]]

st.sidebar.title("Stochastic Bellhop UQ")
st.sidebar.markdown("---")

scen  = st.sidebar.selectbox("Scenario",     SCENARIOS, index=0)
param = st.sidebar.selectbox(
    "Parameter", PARAMS, index=1,
    format_func=lambda p: PARAM_LBLS[p],
)
dist  = st.sidebar.selectbox("Distribution", DISTS, index=2)

default_fom = get_fom(cfg, scen)
fom_val = st.sidebar.number_input(
    "FOM (dB)", value=float(default_fom), step=1.0, min_value=0.0, max_value=200.0
)

st.sidebar.markdown("---")
st.sidebar.subheader("Detection overlay")
thresh = st.sidebar.slider("P(detect) threshold", 0.50, 0.99, 0.70, 0.01)
method = st.sidebar.radio(
    "Overlay method", ["MC", "PCE", "Delta LN3", "None"], index=0
)

# ── Load MC (always required) ─────────────────────────────────────────────────
mc = load_mc(scen, dist, param)

if mc is None:
    st.warning(
        f"No MC results found for **{scen} / {dist} / {param}**. "
        "Run `run_MC.m` first."
    )
    st.stop()

r_km = mc["r_km"]
z_m  = mc["z_m"]
EX   = mc["EX"]         # [Nz × Nr]
Var  = mc["Var"]
P_mc = mc["P_detect"]

# Gaussian P(detect) recomputed at chosen FOM (allows FOM slider without rerunning MATLAB)
EX_safe  = np.where(np.isfinite(EX),  EX,  np.nanmean(EX))
Var_safe = np.where(np.isfinite(Var) & (Var > 0), Var, 1e-6)
P_gauss  = norm_cdf((fom_val - EX_safe) / np.sqrt(Var_safe))

# ── Choose overlay map ────────────────────────────────────────────────────────
if method == "MC":
    P_overlay = P_gauss if abs(fom_val - default_fom) > 0.5 else P_mc
elif method == "PCE":
    pce = load_pce(scen, dist)
    pi  = PARAMS.index(param)
    P_overlay = pce["Pdetect"][pi] if pce is not None else None
elif method == "Delta LN3":
    dlt = load_delta(scen, dist, param)
    P_overlay = dlt["DeltaLN3_prob"] if dlt is not None else None
else:
    P_overlay = None

# Guard: shape mismatch (e.g. PCE/Delta ran at different resolution)
if P_overlay is not None and P_overlay.shape != EX.shape:
    st.sidebar.warning(f"Overlay shape {P_overlay.shape} ≠ EX shape {EX.shape}. Overlay hidden.")
    P_overlay = None


# ══════════════════════════════════════════════════════════════════════════════
# Page title
# ══════════════════════════════════════════════════════════════════════════════
st.title(
    f"TL Map — {scen} | {PARAM_LBLS[param]} | {dist.replace('_', ' ')} | FOM={fom_val:.0f} dB"
)

col_map, col_qtip = st.columns([3, 2])

# ══════════════════════════════════════════════════════════════════════════════
# Left column — heatmap
# ══════════════════════════════════════════════════════════════════════════════
with col_map:
    fig = go.Figure()

    # EX(TL) heatmap
    vmin = float(np.nanpercentile(EX, 2))
    vmax = float(np.nanpercentile(EX, 98))
    fig.add_trace(go.Heatmap(
        x=r_km.tolist(), y=z_m.tolist(), z=EX.tolist(),
        colorscale="Jet",
        zmin=vmin, zmax=vmax,
        colorbar=dict(title="EX(TL) dB", thickness=14, len=0.85),
        name="EX(TL)",
        hovertemplate="r=%.2f km  z=%.1f m  EX=%.1f dB<extra></extra>",
    ))

    # Detection overlay — semi-transparent mask
    if P_overlay is not None:
        det_mask = np.where(P_overlay >= thresh, 1.0, np.nan)
        fig.add_trace(go.Heatmap(
            x=r_km.tolist(), y=z_m.tolist(), z=det_mask.tolist(),
            colorscale=[[0, "rgba(0,220,80,0.0)"], [1, "rgba(0,220,80,0.35)"]],
            showscale=False,
            name=f"Detect ≥ {thresh:.2f}  [{method}]",
            hovertemplate="P(detect)≥{:.2f}<extra></extra>".format(thresh),
        ))
        # Green contour at threshold
        fig.add_trace(go.Contour(
            x=r_km.tolist(), y=z_m.tolist(), z=P_overlay.tolist(),
            contours=dict(start=thresh, end=thresh, size=0, coloring="lines"),
            line=dict(color="lime", width=2),
            showscale=False,
            name="",
            hoverinfo="skip",
        ))

    fig.update_layout(
        xaxis_title="Range (km)",
        yaxis_title="Depth (m)",
        yaxis=dict(autorange="reversed"),
        height=480,
        margin=dict(l=60, r=10, t=30, b=50),
    )
    st.plotly_chart(fig, use_container_width=True, key="main_map")

    # ── Q-TIP coordinate sliders (below map) ──────────────────────────────────
    r_min, r_max = float(r_km.min()), float(r_km.max())
    z_min, z_max = float(z_m.min()),  float(z_m.max())

    c1, c2 = st.columns(2)
    with c1:
        sel_r = st.slider("Range (km)", r_min, r_max,
                          float(np.median(r_km)), step=float((r_max - r_min) / max(len(r_km) - 1, 1)))
    with c2:
        sel_z = st.slider("Depth (m)", z_min, z_max,
                          float(np.median(z_m)),  step=float((z_max - z_min) / max(len(z_m)  - 1, 1)))

# Nearest grid indices
ri = int(np.argmin(np.abs(r_km - sel_r)))
zi = int(np.argmin(np.abs(z_m  - sel_z)))


# ══════════════════════════════════════════════════════════════════════════════
# Right column — Q-TIP panel
# ══════════════════════════════════════════════════════════════════════════════
with col_qtip:
    st.subheader("Q-TIP — Point Inspector")
    st.markdown(f"**r = {r_km[ri]:.2f} km,  z = {z_m[zi]:.1f} m**")

    # ── Stats table ───────────────────────────────────────────────────────────
    rows = []

    # MC row
    ex_px   = float(EX[zi, ri])
    var_px  = float(Var[zi, ri])
    pdet_mc = float(P_mc[zi, ri])
    pdet_mc_fom = float(P_gauss[zi, ri])
    rows.append({
        "Method":      "MC",
        "EX (dB)":     f"{ex_px:.2f}",
        "Var (dB²)":   f"{var_px:.3f}",
        f"P(det|FOM={fom_val:.0f})": f"{pdet_mc_fom:.3f}",
    })

    # Delta rows
    dlt = load_delta(scen, dist, param)
    if dlt is not None:
        def _px(arr, zi, ri):
            a = np.atleast_2d(arr)
            if a.shape[0] > zi and a.shape[1] > ri:
                return float(a[zi, ri])
            return float(np.nan)

        v1 = _px(dlt["Vs_1st"], zi, ri)
        v2 = _px(dlt["Var_TL"], zi, ri)
        pln3 = _px(dlt["DeltaLN3_prob"], zi, ri)
        pd1  = float(norm_cdf((fom_val - ex_px) / max(v1**0.5, 1e-5))) if np.isfinite(v1) else np.nan
        pd2  = float(norm_cdf((fom_val - ex_px) / max(v2**0.5, 1e-5))) if np.isfinite(v2) else np.nan
        rows.append({
            "Method":      "Δ 1st order",
            "EX (dB)":     "—",
            "Var (dB²)":   f"{v1:.3f}" if np.isfinite(v1) else "n/a",
            f"P(det|FOM={fom_val:.0f})": f"{pd1:.3f}" if np.isfinite(pd1) else "n/a",
        })
        rows.append({
            "Method":      "Δ 2nd order",
            "EX (dB)":     "—",
            "Var (dB²)":   f"{v2:.3f}" if np.isfinite(v2) else "n/a",
            f"P(det|FOM={fom_val:.0f})": f"{pd2:.3f}" if np.isfinite(pd2) else "n/a",
        })
        rows.append({
            "Method":      "Δ LN3",
            "EX (dB)":     "—",
            "Var (dB²)":   "—",
            f"P(det|FOM={fom_val:.0f})": f"{pln3:.3f}" if np.isfinite(pln3) else "n/a",
        })

    # PCE row
    pce = load_pce(scen, dist)
    if pce is not None:
        pi     = PARAMS.index(param)
        kstar  = int(pce["Kstar"][pi])
        pd_arr = pce["Pdetect"][pi]
        pd_pce = _px(pd_arr, zi, ri) if dlt is not None else float(np.atleast_2d(pd_arr)[zi, ri])
        rows.append({
            "Method":      f"PCE (K*={kstar})",
            "EX (dB)":     "—",
            "Var (dB²)":   "—",
            f"P(det|FOM={fom_val:.0f})": f"{pd_pce:.3f}" if np.isfinite(pd_pce) else "n/a",
        })

    st.dataframe(
        pd.DataFrame(rows).set_index("Method"),
        use_container_width=True,
    )

    # ── Histogram ─────────────────────────────────────────────────────────────
    st.markdown(f"**MC TL distribution  (Normal 10 pct base samples)**")
    TL_cube = load_tl_cache(scen, param)

    if TL_cube is not None and TL_cube.ndim == 3:
        nz, nr, _ = TL_cube.shape
        if zi < nz and ri < nr:
            tl_px = TL_cube[zi, ri, :]
            fig_h = go.Figure()
            fig_h.add_trace(go.Histogram(
                x=tl_px.tolist(),
                nbinsx=30,
                histnorm="probability density",
                marker_color="steelblue",
                opacity=0.8,
                name="MC",
            ))
            fig_h.add_vline(
                x=fom_val,
                line_color="red", line_dash="dash",
                annotation_text=f"FOM={fom_val:.0f}",
                annotation_position="top right",
            )
            fig_h.update_layout(
                xaxis_title="TL (dB)",
                yaxis_title="Density",
                height=240,
                margin=dict(l=50, r=10, t=20, b=40),
                showlegend=False,
            )
            st.plotly_chart(fig_h, use_container_width=True, key="hist")
        else:
            st.caption("Pixel index out of bounds for TL cache.")
    else:
        st.caption(
            "TL cache not found. Expected at "
            f"`Cache/{scen}/Normal_10pct/TL_{param}.mat`."
        )

    # ── Variance map (EX + Var side by side) ──────────────────────────────────
    st.markdown("**Variance map — MC**")
    fig_var = go.Figure(go.Heatmap(
        x=r_km.tolist(), y=z_m.tolist(), z=Var.tolist(),
        colorscale="Plasma",
        colorbar=dict(title="Var (dB²)", thickness=12, len=0.8),
        hovertemplate="r=%.2f km  z=%.1f m  Var=%.3f<extra></extra>",
    ))
    fig_var.add_trace(go.Scatter(
        x=[r_km[ri]], y=[z_m[zi]],
        mode="markers",
        marker=dict(color="white", size=10, symbol="cross"),
        showlegend=False,
        hoverinfo="skip",
    ))
    fig_var.update_layout(
        xaxis_title="Range (km)",
        yaxis_title="Depth (m)",
        yaxis=dict(autorange="reversed"),
        height=260,
        margin=dict(l=60, r=10, t=20, b=40),
    )
    st.plotly_chart(fig_var, use_container_width=True, key="var_map")
