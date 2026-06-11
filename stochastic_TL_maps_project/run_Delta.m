function run_Delta(sc_target, dist_target)
%% run_Delta — Delta Method UQ across all scenarios and distributions.
%  Call with no args to run all combos, or run_Delta('baseline','Normal_10pct')
%  to process one combo only.
%
% MATHEMATICAL MODEL
% ──────────────────
%   TL(θ) ≈ TL(θ₀) + J'(θ − θ₀)          [1st-order Taylor]
%   E[TL]   ≈ TL(θ₀)
%   Var[TL] ≈ Σᵢ∈S  Jᵢ² σᵢ²               [independent inputs]
%
% Jacobians via central FD  →  computeJacobian.m
% SVP noise                 →  makeSVPNoise.m
%
% Outputs per scenario × distribution:
%   Methods/Delta/<scenario>/<dist>/results/delta_jacobians.mat
%   Methods/Delta/<scenario>/<dist>/results/delta_<subset>.mat  × 7
%   Methods/Delta/<scenario>/<dist>/figures/  (14 PNGs)
%     ├── cheb_<subset>.png         — 4-color categorical Chebyshev map × 7
%     └── smooth_<subset>.png       — non-smooth | LPF | Gaussian × 7

if nargin < 2, sc_target = ''; dist_target = ''; end
close all;
warning('off', 'MATLAB:unknownObjectIEEE');
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:rankDeficientMatrix');
try
    cd(fileparts(mfilename('fullpath')));
catch
end

addpath(genpath('Shared_Utils'));
addpath(genpath('Bellhop'));

cfg = loadConfig();
[subset_names, subset_labels, active] = subsetDefs();

FOM   = cfg.nominal.FOM_dB;
freq0 = cfg.nominal.freq_Hz;
zS0   = cfg.nominal.zS_m;
geo   = cfg.nominal.geo(:)';

h_freq = cfg.jacobian_steps.freq_Hz;
h_zS   = cfg.jacobian_steps.zS_m;
h_svp  = cfg.jacobian_steps.svp_C;

THRESHOLDS = cfg.thresholds(:)';

%% ── LOOP: scenarios × distributions ───────────────────────────────────────
for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);

    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);
        if ~isempty(sc_target) && ~(strcmp(sc.name,sc_target) && strcmp(dist.name,dist_target))
            continue;
        end
        B    = computeVarianceBounds(cfg, dist);

        fprintf('\n════════════════════════════════════════\n');
        fprintf('  Delta | %s | %s\n', sc.name, dist.name);
        fprintf('  σ²_freq=%.2f  σ²_zS=%.6f  σ²_svp=%.6f\n', ...
                B.var_freq, B.var_zS, B.var_svp);
        fprintf('════════════════════════════════════════\n');

        res_dir = fullfile('Methods','Delta', sc.name, dist.name, 'results');
        fig_dir = fullfile('Methods','Delta', sc.name, dist.name, 'figures');
        cache_dir = fullfile('Cache', sc.name, 'bellhop_raw');
        for d = {res_dir, fig_dir, cache_dir}
            if ~exist(d{1},'dir'), mkdir(d{1}); end
        end

        jac_file = fullfile(res_dir, 'delta_jacobians.mat');

        % Bathymetry for this scenario (for overlaying on figures)
        bathy_m = bathymetryMaker(sc.bathy_type, sc.maxR_m);

        %% ── JACOBIANS (only one set per scenario; shared across distributions) ──
        shared_jac = fullfile('Methods','Delta', sc.name, 'shared_jacobians.mat');

        need_compute = ~isfile(shared_jac);
        if ~need_compute
            J = load(shared_jac);
            need_compute = ~isfield(J, 'H_freq');  % recompute if Hessians missing
        end

        if ~need_compute
            fprintf('  [Loading shared Jacobians+Hessians for %s]\n', sc.name);
            TL_nom = J.TL_nom;  r_grid = J.r_grid;  z_grid = J.z_grid;
            r_km = J.r_km;      z_m    = J.z_m;
            J_freq = J.J_freq;  J_zS = J.J_zS;  J_svp = J.J_svp;
            H_freq = J.H_freq;  H_zS = J.H_zS;  H_svp = J.H_svp;
        else
            fprintf('  [Computing Jacobians + diagonal Hessians (2nd-order Delta)]\n');
            max_depth = sc.maxDepth_m;

            sim_nom = {freq0, sc.maxR_m, zS0, 0, 0, 'summer', sc.bathy_type, geo, FOM};

            fprintf('    [1/7] Nominal...\n');
            [TL_nom, r_grid, z_grid] = bellhopCached(sim_nom, cache_dir);
            r_km = r_grid/1000;  z_m = z_grid;

            fprintf('    [2/7] freq+h...\n');
            sim = sim_nom; sim{1} = freq0 + h_freq;
            [TL_fUp,~,~] = bellhopCached(sim, cache_dir);

            fprintf('    [3/7] freq-h...\n');
            sim{1} = freq0 - h_freq;
            [TL_fDn,~,~] = bellhopCached(sim, cache_dir);

            fprintf('    [4/7] zS+h...\n');
            sim = sim_nom; sim{3} = zS0 + h_zS;
            [TL_zUp,~,~] = bellhopCached(sim, cache_dir);

            fprintf('    [5/7] zS-h...\n');
            sim{3} = zS0 - h_zS;
            [TL_zDn,~,~] = bellhopCached(sim, cache_dir);

            fprintf('    [6/7] svp+h...\n');
            sim = sim_nom; sim{6} = 'custom';
            [TL_sUp,~,~] = bellhopCached(sim, cache_dir, ...
                'CustomSVP', makeSVPNoise(+h_svp, max_depth, cfg));

            fprintf('    [7/7] svp-h...\n');
            [TL_sDn,~,~] = bellhopCached(sim, cache_dir, ...
                'CustomSVP', makeSVPNoise(-h_svp, max_depth, cfg));

            % 1st-order Jacobians (central FD)
            J_freq = computeJacobian(TL_fUp, TL_fDn, h_freq);
            J_zS   = computeJacobian(TL_zUp, TL_zDn, h_zS);
            J_svp  = computeJacobian(TL_sUp, TL_sDn, h_svp);

            % 2nd-order diagonal Hessians — FREE (same ±h runs, no new Bellhop calls)
            % H_ii = (TL(θ+h) - 2·TL(θ) + TL(θ-h)) / h²
            H_freq = computeHessianDiag(TL_fUp, TL_nom, TL_fDn, h_freq);
            H_zS   = computeHessianDiag(TL_zUp, TL_nom, TL_zDn, h_zS);
            H_svp  = computeHessianDiag(TL_sUp, TL_nom, TL_sDn, h_svp);

            save(shared_jac, 'TL_nom','J_freq','J_zS','J_svp', ...
                 'H_freq','H_zS','H_svp', ...
                 'r_grid','z_grid','r_km','z_m', '-v7.3');
            fprintf('  Jacobians + Hessians saved to shared cache.\n');
        end

        % 1st-order variance maps per parameter
        V_freq = J_freq.^2 * B.var_freq;
        V_zS   = J_zS.^2   * B.var_zS;
        V_svp  = J_svp.^2  * B.var_svp;

        % 2nd-order Hessian variance coefficient — depends on distribution kurtosis:
        %   Normal:  κᵢ = σᵢ⁴ / 2   (4th central moment = 3σ⁴ → excess = 2σ⁴ → coeff = ½)
        %   Uniform: κᵢ = σᵢ⁴ / 5   (4th central moment = 9σ⁴/5 → excess = 4σ⁴/5 → coeff = ⅕)
        if strcmp(dist.type, 'normal')
            kappa_freq = B.var_freq^2 / 2;
            kappa_zS   = B.var_zS^2   / 2;
            kappa_svp  = B.var_svp^2  / 2;
        else  % uniform
            kappa_freq = B.var_freq^2 / 5;
            kappa_zS   = B.var_zS^2   / 5;
            kappa_svp  = B.var_svp^2  / 5;
        end

        save(jac_file, 'TL_nom','J_freq','J_zS','J_svp', ...
             'H_freq','H_zS','H_svp', ...
             'V_freq','V_zS','V_svp', 'r_grid','z_grid','r_km','z_m', ...
             'FOM','freq0','zS0', 'B', '-v7.3');

        %% ── 7 SUBSETS — 2nd-order Delta Method ────────────────────────────────
        %
        % E[TL]   ≈ TL₀ + ½ Σᵢ Hᵢᵢ σᵢ²           (bias correction)
        % Var[TL] ≈ Σᵢ Jᵢ² σᵢ² + Σᵢ κᵢ Hᵢᵢ²       (curvature-corrected variance)
        %
        Var_delta    = cell(1,7);
        Cheb_lb      = cell(1,7);
        DeltaLN3_map = cell(1,7);
        dE_maps      = cell(1,7);

        for s = 1:7
            fl = active{s};

            % 1st-order variance
            Vs_1st = fl(1)*V_freq + fl(2)*V_zS + fl(3)*V_svp;

            % 2nd-order mean bias correction: ½ Σᵢ Hᵢᵢ σᵢ²
            dE = fl(1)*0.5*H_freq*B.var_freq + ...
                 fl(2)*0.5*H_zS  *B.var_zS   + ...
                 fl(3)*0.5*H_svp *B.var_svp;
            TL_expected = TL_nom + dE;

            % 2nd-order variance: adds κᵢ Hᵢᵢ² per active parameter
            Vs = Vs_1st + fl(1)*H_freq.^2*kappa_freq + ...
                          fl(2)*H_zS.^2  *kappa_zS   + ...
                          fl(3)*H_svp.^2 *kappa_svp;

            Var_delta{s}    = Vs;
            Cheb_lb{s}      = chebyshevBound(TL_expected, Vs, FOM);
            DeltaLN3_map{s} = deltaLN3prob(TL_expected, Vs, FOM);
            dE_maps{s}      = dE;

            sfile = fullfile(res_dir, sprintf('delta_%s.mat', subset_names{s}));
            Var_TL        = Vs;
            Cheb_lb_s     = Cheb_lb{s};
            DeltaLN3_prob = DeltaLN3_map{s};
            save(sfile, 'TL_expected','Var_TL','Cheb_lb_s', ...
                 'Vs_1st','dE','DeltaLN3_prob', ...
                 'r_grid','z_grid','r_km','z_m','FOM');
        end
        fprintf('  All 7 subset results saved.\n');

        %% ── FIGURES ───────────────────────────────────────────────────────────

        % Gaussian blur params (anisotropic, physically motivated)
        dz = z_m(2) - z_m(1);             % m/pixel (depth)
        dr = (r_km(2) - r_km(1)) * 1000;  % m/pixel (range)
        SZ_M  = 1.0;    % sigma_z in metres  (light depth smoothing)
        SR_M  = 100;    % sigma_r in metres  (light range smoothing)
        sz_px = SZ_M / dz;
        sr_px = SR_M / dr;
        % LPF: anisotropic cutoffs derived from physical feature size
        % cutoff = 2 * pixel_spacing / min_feature_size_metres
        MIN_FEAT_Z = 20;    % filter depth features smaller than 20 m
        MIN_FEAT_R = 500;   % filter range features smaller than 500 m
        lpf_cut_z  = min(0.90, 2 * dz / MIN_FEAT_Z);
        lpf_cut_r  = min(0.90, 2 * dr / MIN_FEAT_R);

        for s = 1:7
            lbl = subset_labels{s};
            C   = Cheb_lb{s};

            %% Fig A: 4-color Chebyshev probability bound map
            figA = figure('Position',[50 50 700 480]);
            shadowCategoryMap(gca, r_km, z_m, C, ...
                sprintf('Delta Chebyshev | %s | %s | %s', sc.name, dist.name, lbl), ...
                THRESHOLDS);
            overlayBathymetry(gca, bathy_m, sc.maxDepth_m);
            saveFigPNG(figA, fullfile(fig_dir, sprintf('cheb_%s', subset_names{s})));
            drawnow; close(figA);

            %% Fig A2: Delta-LN3 probability map (moment-matched, no MC)
            figA2  = figure('Position',[50 50 700 480]);
            shadowCategoryMap(gca, r_km, z_m, DeltaLN3_map{s}, ...
                sprintf('Delta-LN3 | %s | %s | %s', sc.name, dist.name, lbl), ...
                THRESHOLDS);
            overlayBathymetry(gca, bathy_m, sc.maxDepth_m);
            saveFigPNG(figA2, fullfile(fig_dir, sprintf('cheb_ln3_%s', subset_names{s})));
            drawnow; close(figA2);

            %% Fig B: Smoothing comparison (3-panel) + individual smoothed PNGs
            C_lpf  = lpf2d(C, lpf_cut_z, lpf_cut_r);
            C_gaus = gaussBlur2d(C, sz_px, sr_px);

            lpf_lbl  = sprintf('LPF (z>%.0fm, r>%.0fm)', MIN_FEAT_Z, MIN_FEAT_R);
            gaus_lbl = sprintf('Gaussian (sz=%.0fm, sr=%.0fm)', SZ_M, SR_M);

            % 3-panel comparison: Non-smooth | LPF | Gaussian
            figB = figure('Position',[50 50 1600 480]);
            ax1 = subplot(1,3,1);
            shadowCategoryMap(ax1, r_km, z_m, C,      'Non-smoothed', THRESHOLDS);
            overlayBathymetry(ax1, bathy_m, sc.maxDepth_m);
            ax2 = subplot(1,3,2);
            shadowCategoryMap(ax2, r_km, z_m, C_lpf,  lpf_lbl,        THRESHOLDS);
            overlayBathymetry(ax2, bathy_m, sc.maxDepth_m);
            ax3 = subplot(1,3,3);
            shadowCategoryMap(ax3, r_km, z_m, C_gaus, gaus_lbl,       THRESHOLDS);
            overlayBathymetry(ax3, bathy_m, sc.maxDepth_m);
            sgtitle(sprintf('Delta Smoothing | %s | %s | %s | FOM=%ddB', ...
                            sc.name, dist.name, lbl, FOM));
            saveFigPNG(figB, fullfile(fig_dir, sprintf('smooth_%s', subset_names{s})));
            drawnow; close(figB);

            % Individual smoothed PNGs (LPF)
            figL = figure('Position',[50 50 700 480]);
            shadowCategoryMap(gca, r_km, z_m, C_lpf, ...
                sprintf('Delta %s | %s | %s | %s', lpf_lbl, sc.name, dist.name, lbl), ...
                THRESHOLDS);
            overlayBathymetry(gca, bathy_m, sc.maxDepth_m);
            saveFigPNG(figL, fullfile(fig_dir, sprintf('smooth_lpf_%s', subset_names{s})));
            drawnow; close(figL);

            % Individual smoothed PNGs (Gaussian)
            figG = figure('Position',[50 50 700 480]);
            shadowCategoryMap(gca, r_km, z_m, C_gaus, ...
                sprintf('Delta %s | %s | %s | %s', gaus_lbl, sc.name, dist.name, lbl), ...
                THRESHOLDS);
            overlayBathymetry(gca, bathy_m, sc.maxDepth_m);
            saveFigPNG(figG, fullfile(fig_dir, sprintf('smooth_gauss_%s', subset_names{s})));
            drawnow; close(figG);
        end

        %% Fig: 2nd-order bias map (all 7 subsets in one figure)
        % Shows dE = ½ Σ Hᵢᵢ σᵢ² — the mean TL shift from curvature
        dE_abs_max = max(cellfun(@(m) max(abs(m(:))), dE_maps));
        dE_abs_max = max(dE_abs_max, 0.1);  % avoid zero clim

        fw = max(1400, 200*7);
        fig_bias = figure('Position',[50 50 fw 420]);
        for s = 1:7
            ax = subplot(1,7,s);
            imagesc(r_km, z_m, dE_maps{s});
            colormap(ax, redblue(256)); colorbar(ax);
            clim([-dE_abs_max, dE_abs_max]);
            set(ax,'YDir','reverse');
            xlabel(ax,'Range (km)'); if s==1, ylabel(ax,'Depth (m)'); end
            title(ax, subset_labels{s}, 'FontSize',8,'Interpreter','none');
            overlayBathymetry(ax, bathy_m, sc.maxDepth_m);
        end
        sgtitle(sprintf('2nd-order Delta bias  ½ΣHᵢᵢσᵢ²  |  %s | %s  (dB)', ...
                        sc.name, dist.name));
        saveFigPNG(fig_bias, fullfile(fig_dir,'bias_2nd_order'));
        drawnow; close(fig_bias);

        fprintf('  Figures saved: %s\n', fig_dir);
    end
end

fprintf('\n=== run_Delta complete. ===\n');
end
