function compute_delta_var_tables()
%% compute_delta_var_tables  Var-deviation summary tables (MC vs Delta).
%
% Loads existing results only (no Bellhop, no recomputation):
%   Methods/MC/<sc>/<dist>/results/MC_<param>.mat      → MC_Var
%   Methods/Delta/<sc>/<dist>/results/delta_<param>.mat → Vs_1st (1st-order
%       Delta variance J²σ²) and Var_TL (2nd-order total = Vs_1st + κH²)
%   (field names per pipeline/run_delta.m, single-parameter subsets)
%
% TABLE A — mean absolute deviation across ALL map pixels between MC's Var
%   and Delta's 1st-order Var, mean(|MC_Var - Vs_1st|) in dB², averaged
%   over the 3 params (freq, zS, svp).  Rows = scenarios, cols = dists.
%
% TABLE B — average 2nd-order Hessian correction magnitude relative to the
%   1st-order variance: 100 * mean(Var_TL - Vs_1st) / mean(Vs_1st) in %,
%   for all (param × distribution) combos, averaged over the 4 scenarios.
%   (Var_TL - Vs_1st = κᵢ Hᵢᵢ² exactly, per run_delta.m.)
%
% Prints both as ready-to-paste LaTeX tabular snippets (3 sig. figs),
% delimited by %%% TABLE A %%% / %%% TABLE B %%%, and saves the raw
% matrices to Methods/Delta/var_deviation_tables.mat.

ROOT = fileparts(fileparts(mfilename('fullpath')));
try; cd(ROOT); catch; end
addpath(genpath(fullfile(ROOT,'core')));

cfg    = loadConfig();
SCEN   = {cfg.scenarios.name};          % 4
DISTS  = {cfg.distributions.name};      % 3
PARAMS = {'freq','zS','svp'};
nS = numel(SCEN);  nD = numel(DISTS);  nP = numel(PARAMS);

% mad_A(si,di,pi)  = mean(|MC_Var - Vs_1st|) over pixels (dB²)
% rel_B(pi,di,si)  = 100*mean(corr2nd)/mean(Vs_1st) (%)
mad_A = NaN(nS, nD, nP);
rel_B = NaN(nP, nD, nS);

for si = 1:nS
    for di = 1:nD
        for pi = 1:nP
            mc_f = fullfile('Methods','MC', SCEN{si}, DISTS{di}, 'results', ...
                            sprintf('MC_%s.mat', PARAMS{pi}));
            dl_f = fullfile('Methods','Delta', SCEN{si}, DISTS{di}, 'results', ...
                            sprintf('delta_%s.mat', PARAMS{pi}));
            if ~isfile(mc_f), fprintf('[MISS] %s\n', mc_f); continue; end
            if ~isfile(dl_f), fprintf('[MISS] %s\n', dl_f); continue; end

            mc = load(mc_f, 'MC_Var');
            dl = load(dl_f, 'Vs_1st', 'Var_TL');

            V_mc = double(mc.MC_Var(:));
            V1   = double(dl.Vs_1st(:));
            V2   = double(dl.Var_TL(:));
            ok   = isfinite(V_mc) & isfinite(V1) & isfinite(V2);

            mad_A(si,di,pi) = mean(abs(V_mc(ok) - V1(ok)));
            corr2 = V2(ok) - V1(ok);              % = kappa * H^2 exactly
            rel_B(pi,di,si) = 100 * mean(corr2) / max(mean(V1(ok)), eps);
        end
    end
end

A = mean(mad_A, 3, 'omitnan');   % [nS x nD] averaged over params
B = mean(rel_B, 3, 'omitnan');   % [nP x nD] averaged over scenarios

out_mat = fullfile('Methods','Delta','var_deviation_tables.mat');
save(out_mat, 'mad_A','rel_B','A','B','SCEN','DISTS','PARAMS');
fprintf('Saved raw matrices: %s\n\n', out_mat);

% Literal LaTeX labels (printed via %s — no printf re-escaping)
dist_lbl = {'1\%','5\%','10\%'};
scen_lbl = strrep(SCEN, '_', '\_');

%% ── TABLE A ──────────────────────────────────────────────────────────────────
fprintf('%%%%%% TABLE A %%%%%%\n');
fprintf('%% Mean absolute deviation |Var_MC - Var_Delta(1st)| over all map pixels,\n');
fprintf('%% averaged over the 3 params (freq, zS, svp).  Units: dB^2.\n');
fprintf('\\begin{tabular}{lccc}\n\\hline\n');
fprintf('Scenario & Normal %s & Normal %s & Normal %s \\\\\n', dist_lbl{:});
fprintf('\\hline\n');
for si = 1:nS
    fprintf('%s', scen_lbl{si});
    for di = 1:nD
        fprintf(' & %s', sig3(A(si,di)));
    end
    fprintf(' \\\\\n');
end
fprintf('\\hline\n\\end{tabular}\n');
fprintf('%%%%%% END TABLE A %%%%%%\n\n');

%% ── TABLE B ──────────────────────────────────────────────────────────────────
fprintf('%%%%%% TABLE B %%%%%%\n');
fprintf('%% Mean 2nd-order Hessian correction kappa*H^2 relative to 1st-order\n');
fprintf('%% variance J^2*sigma^2, 100*mean(corr)/mean(Var1) in %%, averaged over\n');
fprintf('%% all pixels and all 4 scenarios.\n');
fprintf('\\begin{tabular}{lccc}\n\\hline\n');
fprintf('Parameter & Normal %s & Normal %s & Normal %s \\\\\n', dist_lbl{:});
fprintf('\\hline\n');
param_tex = {'freq','$z_S$','SVP'};
for pi = 1:nP
    fprintf('%s', param_tex{pi});
    for di = 1:nD
        fprintf(' & %s\\%%', sig3(B(pi,di)));
    end
    fprintf(' \\\\\n');
end
fprintf('\\hline\n\\end{tabular}\n');
fprintf('%%%%%% END TABLE B %%%%%%\n');
end

function s = sig3(v)
% Format to 3 significant figures.
    if ~isfinite(v), s = '--'; return; end
    s = sprintf('%.3g', v);
end
