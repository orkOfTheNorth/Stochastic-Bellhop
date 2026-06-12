function run_Delta(sc_target, dist_target)
%% run_Delta — Delta Method UQ across all scenarios and distributions.
%  Call with no args to run all combos, or run_Delta('deep_water','Normal_10pct')
%  to process one combo only.
%
% MATHEMATICAL MODEL
% ──────────────────
%   TL(θ) ≈ TL(θ₀) + J'(θ − θ₀)           [1st-order Taylor]
%   Var[TL] ≈ Σᵢ∈S  Jᵢ² σᵢ²                [independent inputs]
%   Var[TL] ≈ Σᵢ    Jᵢ² σᵢ² + κᵢ Hᵢᵢ²     [2nd-order, diagonal Hessian]
%
% GH ORDER SWEEP (new — orders 1-10)
% ────────────────────────────────────
%   For K = 1..10: variance estimated via K-point Gauss-Hermite quadrature
%   per parameter (σ-scaled nodes, progressively sorted by weight).
%   Total extra Bellhop runs: 10 × 3 params = 30 per scenario × distribution.
%   All go through bellhopCached (cached on disk).
%
% Outputs per scenario × distribution:
%   Methods/Delta/<scenario>/<dist>/results/delta_jacobians.mat
%   Methods/Delta/<scenario>/<dist>/results/delta_<subset>.mat  × 7
%   Methods/Delta/<scenario>/<dist>/results/delta_GH_sweep.mat
%   Methods/Delta/<scenario>/<dist>/figures/  (≥20 PNGs)

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

freq0 = cfg.nominal.freq_Hz;
zS0   = cfg.nominal.zS_m;
geo   = cfg.nominal.geo(:)';

h_freq = cfg.jacobian_steps.freq_Hz;
h_zS   = cfg.jacobian_steps.zS_m;
h_svp  = cfg.jacobian_steps.svp_C;

THRESHOLDS    = cfg.thresholds(:)';
GH_MAX_ORDER  = 10;

%% ── LOOP: scenarios × distributions ────────────────────────────────────────
for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);

    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);
        if ~isempty(sc_target) && ~(strcmp(sc.name,sc_target) && strcmp(dist.name,dist_target))
            continue;
        end

        FOM = getFOM(cfg, sc.name);
        B   = computeVarianceBounds(cfg, dist);

        fprintf('\n════════════════════════════════════════\n');
        fprintf('  Delta | %s | %s | FOM=%ddB\n', sc.name, dist.name, FOM);
        fprintf('  σ²_freq=%.2f  σ²_zS=%.6f  σ²_svp=%.6f\n', ...
                B.var_freq, B.var_zS, B.var_svp);
        fprintf('════════════════════════════════════════\n');

        res_dir   = fullfile('Methods','Delta', sc.name, dist.name, 'results');
        fig_dir   = fullfile('Methods','Delta', sc.name, dist.name, 'figures');
        cache_dir = fullfile('Cache', sc.name, 'bellhop_raw');
        for d = {res_dir, fig_dir, cache_dir}
            if ~exist(d{1},'dir'), mkdir(d{1}); end
        end

        jac_file = fullfile(res_dir, 'delta_jacobians.mat');
        bathy_m  = bathymetryMaker(sc.bathy_type, sc.maxR_m);

        %% ── JACOBIANS (shared across distributions for this scenario) ──────
        shared_jac  = fullfile('Methods','Delta', sc.name, 'shared_jacobians.mat');
        need_compute = ~isfile(shared_jac);
        if ~need_compute
            J            = load(shared_jac);
            need_compute = ~isfield(J, 'H_freq');
        end

        if ~need_compute
            fprintf('  [Loading shared Jacobians+Hessians for %s]\n', sc.name);
            TL_nom = J.TL_nom;  r_grid = J.r_grid;  z_grid = J.z_grid;
            r_km   = J.r_km;    z_m    = J.z_m;
            J_freq = J.J_freq;  J_zS   = J.J_zS;    J_svp  = J.J_svp;
            H_freq = J.H_freq;  H_zS   = J.H_zS;    H_svp  = J.H_svp;
        else
            fprintf('  [Computing Jacobians + diagonal Hessians (2nd-order Delta)]\n');

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
                'CustomSVP', makeSVPNoise(+h_svp, sc.maxDepth_m, cfg));

            fprintf('    [7/7] svp-h...\n');
            [TL_sDn,~,~] = bellhopCached(sim, cache_dir, ...
                'CustomSVP', makeSVPNoise(-h_svp, sc.maxDepth_m, cfg));

            J_freq = computeJacobian(TL_fUp, TL_fDn, h_freq);
            J_zS   = computeJacobian(TL_zUp, TL_zDn, h_zS);
            J_svp  = computeJacobian(TL_sUp, TL_sDn, h_svp);

            H_freq = computeHessianDiag(TL_fUp, TL_nom, TL_fDn, h_freq);
            H_zS   = computeHessianDiag(TL_zUp, TL_nom, TL_zDn, h_zS);
            H_svp  = computeHessianDiag(TL_sUp, TL_nom, TL_sDn, h_svp);

            save(shared_jac, 'TL_nom','J_freq','J_zS','J_svp', ...
                 'H_freq','H_zS','H_svp', ...
                 'r_grid','z_grid','r_km','z_m', '-v7.3');
            fprintf('  Jacobians + Hessians saved to shared cache.\n');
        end

        % sim_nom always available for GH sweep (uses per-loop FOM)
        sim_nom = {freq0, sc.maxR_m, zS0, 0, 0, 'summer', sc.bathy_type, geo, FOM};

        % 1st-order variance maps per parameter
        V_freq = J_freq.^2 * B.var_freq;
        V_zS   = J_zS.^2   * B.var_zS;
        V_svp  = J_svp.^2  * B.var_svp;

        % 2nd-order kappa: depends on distribution type
        %   Normal:  E[ε⁴] = 3σ⁴  →  Var correction = (H/2)² × 2σ⁴ = H²σ⁴/2
        %   Uniform: E[ε⁴] = 9σ⁴/5 →  correction = (H/2)² × 4σ⁴/5 = H²σ⁴/5
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

        %% ── 7 SUBSETS — 2nd-order Delta ─────────────────────────────────────
        %
        % E[TL]   ≈ TL₀ + ½ Σᵢ Hᵢᵢ σᵢ²   (bias correction)
        % Var[TL] ≈ Σᵢ Jᵢ² σᵢ² + Σᵢ κᵢ Hᵢᵢ²
        %
        Var_delta    = cell(1,7);
        Cheb_lb      = cell(1,7);
        DeltaLN3_map = cell(1,7);   % P(detect) = P(TL < FOM)
        dE_maps      = cell(1,7);

        for s = 1:7
            fl = active{s};

            Vs_1st = fl(1)*V_freq + fl(2)*V_zS + fl(3)*V_svp;

            dE = fl(1)*0.5*H_freq*B.var_freq + ...
                 fl(2)*0.5*H_zS  *B.var_zS   + ...
                 fl(3)*0.5*H_svp *B.var_svp;
            TL_expected = TL_nom + dE;

            Vs = Vs_1st + fl(1)*H_freq.^2*kappa_freq + ...
                          fl(2)*H_zS.^2  *kappa_zS   + ...
                          fl(3)*H_svp.^2 *kappa_svp;

            Var_delta{s}    = Vs;
            Cheb_lb{s}      = chebyshevBound(TL_expected, Vs, FOM);
            % Convert P(shadow) → P(detect) = P(TL < FOM)
            DeltaLN3_map{s} = 1 - deltaLN3prob(TL_expected, Vs, FOM);
            dE_maps{s}      = dE;

            sfile         = fullfile(res_dir, sprintf('delta_%s.mat', subset_names{s}));
            Var_TL        = Vs;
            Cheb_lb_s     = Cheb_lb{s};           % P(shadow) lower bound — kept for run_Comparison
            DeltaLN3_prob = DeltaLN3_map{s};      % P(detect) — used by run_Comparison
            save(sfile, 'TL_expected','Var_TL','Cheb_lb_s', ...
                 'Vs_1st','dE','DeltaLN3_prob', ...
                 'r_grid','z_grid','r_km','z_m','FOM');
        end
        fprintf('  All 7 subset results saved.\n');

        %% ── GH ORDER SWEEP (orders 1-10) ────────────────────────────────────
        % For order K: variance estimated via K-point GH quadrature per parameter.
        % Uses GH_MAX_ORDER=10 evaluation points (sorted descending by weight);
        % order K uses the first K of these progressively — no extra runs needed.
        %
        gh_sweep_file = fullfile(res_dir, 'delta_GH_sweep.mat');
        if ~isfile(gh_sweep_file)
            fprintf('  [GH sweep] Orders 1–%d | %s | %s\n', GH_MAX_ORDER, sc.name, dist.name);

            sig_freq = sqrt(B.var_freq);
            sig_zS   = sqrt(B.var_zS);
            sig_svp  = sqrt(B.var_svp);

            [Nz, Nr] = size(TL_nom);

            % Get GH_MAX_ORDER nodes, sort by weight (descend) for progressive use
            [xi_full, w_full] = gaussQuadNodes(GH_MAX_ORDER, 'normal');
            [~, sort_idx]     = sort(w_full, 'descend');
            xi_s = xi_full(sort_idx);   % [GH_MAX_ORDER × 1] sorted nodes
            w_s  = w_full(sort_idx);    % [GH_MAX_ORDER × 1] sorted weights

            % ── Run Bellhop at all GH_MAX_ORDER node locations per parameter ──
            fprintf('    Running Bellhop at %d freq nodes...\n', GH_MAX_ORDER);
            TL_GH_freq = zeros(Nz, Nr, GH_MAX_ORDER);
            for k = 1:GH_MAX_ORDER
                theta_k  = max(freq0 + sig_freq * xi_s(k), 50);
                sim_k    = sim_nom; sim_k{1} = theta_k;
                TL_GH_freq(:,:,k) = bellhopCached(sim_k, cache_dir);
            end

            fprintf('    Running Bellhop at %d zS nodes...\n', GH_MAX_ORDER);
            TL_GH_zS = zeros(Nz, Nr, GH_MAX_ORDER);
            for k = 1:GH_MAX_ORDER
                theta_k = zS0 + sig_zS * xi_s(k);
                theta_k = max(theta_k, 1);
                theta_k = min(theta_k, sc.maxDepth_m * 0.95);
                sim_k   = sim_nom; sim_k{3} = theta_k;
                TL_GH_zS(:,:,k) = bellhopCached(sim_k, cache_dir);
            end

            fprintf('    Running Bellhop at %d SVP nodes...\n', GH_MAX_ORDER);
            TL_GH_svp = zeros(Nz, Nr, GH_MAX_ORDER);
            for k = 1:GH_MAX_ORDER
                pert_k  = sig_svp * xi_s(k);
                sim_k   = sim_nom; sim_k{6} = 'custom';
                TL_GH_svp(:,:,k) = bellhopCached(sim_k, cache_dir, ...
                    'CustomSVP', makeSVPNoise(pert_k, sc.maxDepth_m, cfg));
            end

            % ── Compute Var[TL] at each GH order K (using first K sorted nodes) ──
            Var_GH_freq = zeros(Nz, Nr, GH_MAX_ORDER);
            EX_GH_freq  = zeros(Nz, Nr, GH_MAX_ORDER);
            Var_GH_zS   = zeros(Nz, Nr, GH_MAX_ORDER);
            EX_GH_zS    = zeros(Nz, Nr, GH_MAX_ORDER);
            Var_GH_svp  = zeros(Nz, Nr, GH_MAX_ORDER);
            EX_GH_svp   = zeros(Nz, Nr, GH_MAX_ORDER);

            for K = 1:GH_MAX_ORDER
                w_K  = w_s(1:K) / sum(w_s(1:K));    % normalized sub-weights
                w3   = reshape(w_K, 1, 1, K);

                % Frequency
                EX_K = sum(TL_GH_freq(:,:,1:K) .* w3, 3);
                EX_GH_freq(:,:,K)  = EX_K;
                Var_GH_freq(:,:,K) = sum(w3 .* (TL_GH_freq(:,:,1:K) - EX_K).^2, 3);

                % zS
                EX_K = sum(TL_GH_zS(:,:,1:K) .* w3, 3);
                EX_GH_zS(:,:,K)  = EX_K;
                Var_GH_zS(:,:,K) = sum(w3 .* (TL_GH_zS(:,:,1:K) - EX_K).^2, 3);

                % SVP
                EX_K = sum(TL_GH_svp(:,:,1:K) .* w3, 3);
                EX_GH_svp(:,:,K)  = EX_K;
                Var_GH_svp(:,:,K) = sum(w3 .* (TL_GH_svp(:,:,1:K) - EX_K).^2, 3);
            end

            % Delta 1st- and 2nd-order per-parameter Var (for reference in figures)
            V1_freq = V_freq;
            V1_zS   = V_zS;
            V1_svp  = V_svp;
            V2_freq = V_freq + H_freq.^2 * kappa_freq;
            V2_zS   = V_zS   + H_zS.^2   * kappa_zS;
            V2_svp  = V_svp  + H_svp.^2  * kappa_svp;

            save(gh_sweep_file, ...
                 'Var_GH_freq','Var_GH_zS','Var_GH_svp', ...
                 'EX_GH_freq','EX_GH_zS','EX_GH_svp', ...
                 'TL_GH_freq','TL_GH_zS','TL_GH_svp', ...
                 'V1_freq','V1_zS','V1_svp', ...
                 'V2_freq','V2_zS','V2_svp', ...
                 'xi_s','w_s','sig_freq','sig_zS','sig_svp', ...
                 'r_km','z_m','FOM', '-v7.3');
            fprintf('  GH sweep saved → %s\n', gh_sweep_file);
        else
            fprintf('  [GH sweep] Loaded from cache: %s\n', gh_sweep_file);
            GH = load(gh_sweep_file, 'Var_GH_freq','Var_GH_zS','Var_GH_svp', ...
                      'EX_GH_freq','EX_GH_zS','EX_GH_svp', ...
                      'V1_freq','V1_zS','V1_svp','V2_freq','V2_zS','V2_svp');
            Var_GH_freq = GH.Var_GH_freq;  Var_GH_zS = GH.Var_GH_zS;  Var_GH_svp = GH.Var_GH_svp;
            EX_GH_freq  = GH.EX_GH_freq;   EX_GH_zS  = GH.EX_GH_zS;   EX_GH_svp  = GH.EX_GH_svp;
            V1_freq = GH.V1_freq;  V1_zS = GH.V1_zS;  V1_svp = GH.V1_svp;
            V2_freq = GH.V2_freq;  V2_zS = GH.V2_zS;  V2_svp = GH.V2_svp;
        end

        % Load MC Var for L2 comparison (optional — skip figures if not available)
        mc_base_dir = fullfile('Methods','MC', sc.name, dist.name, 'results');
        mc_ref_file = fullfile(mc_base_dir, 'MC_zS_freq_svp.mat');
        if ~isfile(mc_ref_file)
            mc_ref_file = fullfile(mc_base_dir, 'MC_zS.mat');   % fallback
        end
        have_mc = isfile(mc_ref_file);
        if have_mc
            tmp = load(mc_ref_file, 'MC_Var','MC_EX');
            MC_Var_ref = tmp.MC_Var;
            MC_EX_ref  = tmp.MC_EX;
            clear tmp;
        end

        %% ── FIGURES ─────────────────────────────────────────────────────────

        % Anisotropic Gaussian blur params
        dz = z_m(2) - z_m(1);
        dr = (r_km(2) - r_km(1)) * 1000;
        SZ_M  = 1.0;   SR_M  = 100;
        sz_px = SZ_M / dz;
        sr_px = SR_M / dr;
        MIN_FEAT_Z = 20;    MIN_FEAT_R = 500;
        lpf_cut_z  = min(0.90, 2 * dz / MIN_FEAT_Z);
        lpf_cut_r  = min(0.90, 2 * dr / MIN_FEAT_R);

        for s = 1:7
            lbl = subset_labels{s};
            C   = Cheb_lb{s};          % P(shadow) lower bound
            Pd  = 1 - C;               % upper bound on P(detect)

            %% Fig A: Chebyshev detection bound (inverted to P_detect)
            figA = figure('Position',[50 50 700 480]);
            detectionCategoryMap(gca, r_km, z_m, Pd, ...
                sprintf('Delta Cheb P(detect) UB | %s | %s | %s', sc.name, dist.name, lbl), ...
                THRESHOLDS);
            overlayBathymetry(gca, bathy_m, sc.maxDepth_m);
            saveFigPNG(figA, fullfile(fig_dir, sprintf('cheb_%s', subset_names{s})));
            drawnow; close(figA);

            %% Fig A2: Delta-LN3 P(detect) map
            figA2 = figure('Position',[50 50 700 480]);
            detectionCategoryMap(gca, r_km, z_m, DeltaLN3_map{s}, ...
                sprintf('Delta-LN3 P(detect) | %s | %s | %s', sc.name, dist.name, lbl), ...
                THRESHOLDS);
            overlayBathymetry(gca, bathy_m, sc.maxDepth_m);
            saveFigPNG(figA2, fullfile(fig_dir, sprintf('cheb_ln3_%s', subset_names{s})));
            drawnow; close(figA2);

            %% Fig B: Smoothing comparison (3-panel)
            C_lpf  = lpf2d(Pd, lpf_cut_z, lpf_cut_r);
            C_gaus = gaussBlur2d(Pd, sz_px, sr_px);
            lpf_lbl  = sprintf('LPF (z>%.0fm, r>%.0fm)', MIN_FEAT_Z, MIN_FEAT_R);
            gaus_lbl = sprintf('Gaussian (sz=%.0fm, sr=%.0fm)', SZ_M, SR_M);

            figB = figure('Position',[50 50 1600 480]);
            ax1 = subplot(1,3,1);
            detectionCategoryMap(ax1, r_km, z_m, Pd,     'Non-smoothed',  THRESHOLDS);
            overlayBathymetry(ax1, bathy_m, sc.maxDepth_m);
            ax2 = subplot(1,3,2);
            detectionCategoryMap(ax2, r_km, z_m, C_lpf,  lpf_lbl,         THRESHOLDS);
            overlayBathymetry(ax2, bathy_m, sc.maxDepth_m);
            ax3 = subplot(1,3,3);
            detectionCategoryMap(ax3, r_km, z_m, C_gaus, gaus_lbl,        THRESHOLDS);
            overlayBathymetry(ax3, bathy_m, sc.maxDepth_m);
            sgtitle(sprintf('Delta P(detect) Smoothing | %s | %s | %s | FOM=%ddB', ...
                            sc.name, dist.name, lbl, FOM));
            saveFigPNG(figB, fullfile(fig_dir, sprintf('smooth_%s', subset_names{s})));
            drawnow; close(figB);

            % Individual smoothed maps
            figL = figure('Position',[50 50 700 480]);
            detectionCategoryMap(gca, r_km, z_m, C_lpf, ...
                sprintf('Delta %s | %s | %s | %s', lpf_lbl, sc.name, dist.name, lbl), ...
                THRESHOLDS);
            overlayBathymetry(gca, bathy_m, sc.maxDepth_m);
            saveFigPNG(figL, fullfile(fig_dir, sprintf('smooth_lpf_%s', subset_names{s})));
            drawnow; close(figL);

            figG = figure('Position',[50 50 700 480]);
            detectionCategoryMap(gca, r_km, z_m, C_gaus, ...
                sprintf('Delta %s | %s | %s | %s', gaus_lbl, sc.name, dist.name, lbl), ...
                THRESHOLDS);
            overlayBathymetry(gca, bathy_m, sc.maxDepth_m);
            saveFigPNG(figG, fullfile(fig_dir, sprintf('smooth_gauss_%s', subset_names{s})));
            drawnow; close(figG);
        end

        %% Fig: 2nd-order bias map (all 7 subsets)
        dE_abs_max = max(cellfun(@(m) max(abs(m(:))), dE_maps));
        dE_abs_max = max(dE_abs_max, 0.1);

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

        %% ── GH SWEEP FIGURES ────────────────────────────────────────────────

        param_names = {'freq', 'zS', 'svp'};
        param_lbls  = {'Frequency', 'Source depth z_S', 'SVP shift'};
        Var_GH_all  = {Var_GH_freq, Var_GH_zS, Var_GH_svp};
        V1_all = {V1_freq, V1_zS, V1_svp};
        V2_all = {V2_freq, V2_zS, V2_svp};

        tag = sprintf('%s | %s', sc.name, dist.name);

        for pi = 1:3
            pname    = param_names{pi};
            plbl     = param_lbls{pi};
            VGH_p    = Var_GH_all{pi};    % [Nz × Nr × 10]
            V1_p     = V1_all{pi};
            V2_p     = V2_all{pi};

            %% Fig: 10-panel scatter — GH order 1..10 vs MC Var (or vs GH-10 if MC missing)
            if have_mc
                Vref = MC_Var_ref(:);
                ref_lbl = 'MC Var';
            else
                Vref = VGH_p(:,:,GH_MAX_ORDER);   % GH-10 as reference
                Vref = Vref(:);
                ref_lbl = 'GH order 10';
            end
            n_panels = GH_MAX_ORDER + 2;  % 10 GH orders + Delta-1 + Delta-2
            nc = 4;  nr_fig = ceil(n_panels / nc);
            fig_sc = figure('Position',[20 20 nc*340 nr_fig*300]);

            all_V_series = [V1_p(:), V2_p(:), reshape(VGH_p, [], GH_MAX_ORDER)];
            all_lbls_sc  = [{'Delta-1','Delta-2'}, ...
                             arrayfun(@(k) sprintf('GH K=%d',k), 1:GH_MAX_ORDER, 'UniformOutput',false)];

            L2_orders = zeros(1, n_panels);
            for oi = 1:n_panels
                Vest = all_V_series(:, oi);

                ax = subplot(nr_fig, nc, oi);
                scatter(ax, Vref(:), Vest(:), 2, [0.2 0.4 0.8], 'filled', ...
                        'MarkerFaceAlpha', 0.3);
                hold(ax,'on');
                v_max = max(max(Vref(:)), max(Vest(:)));
                plot(ax, [0 v_max], [0 v_max], 'r-', 'LineWidth',1.2);
                hold(ax,'off');

                % R² (on pixels with Vref > 0)
                vr = Vref(Vref>0);  ve = Vest(Vref>0);
                if numel(vr) > 2
                    ss_tot = sum((vr - mean(vr)).^2);
                    ss_res = sum((vr - ve).^2);
                    R2 = max(1 - ss_res / max(ss_tot, eps), -10);
                else
                    R2 = NaN;
                end
                L2_orders(oi) = sqrt(mean((Vref - Vest).^2));

                title(ax, sprintf('%s  R²=%.3f', all_lbls_sc{oi}, R2), ...
                      'FontSize',8,'Interpreter','none');
                xlabel(ax, ref_lbl, 'FontSize',7);
                ylabel(ax, 'Delta Var', 'FontSize',7);
                axis(ax,'tight');
            end
            sgtitle(sprintf('Delta Var convergence — %s | %s | %s', plbl, tag, ref_lbl), ...
                    'Interpreter','none', 'FontSize',11);
            saveFigPNG(fig_sc, fullfile(fig_dir, sprintf('delta_GH_scatter_%s', pname)));
            drawnow; close(fig_sc);

            %% Fig: L2 error vs GH order (with Delta-1 and Delta-2 as reference lines)
            fig_L2 = figure('Position',[50 50 800 440]);
            orders_x = 1:GH_MAX_ORDER;
            L2_GH = zeros(1, GH_MAX_ORDER);
            for K = 1:GH_MAX_ORDER
                L2_GH(K) = sqrt(mean((Vref - reshape(VGH_p(:,:,K),[],1)).^2));
            end
            L2_d1 = sqrt(mean((Vref - V1_p(:)).^2));
            L2_d2 = sqrt(mean((Vref - V2_p(:)).^2));

            plot(orders_x, L2_GH, 'o-b', 'LineWidth',1.5, 'MarkerSize',7, ...
                 'DisplayName','GH order K'); hold on;
            yline(L2_d1, '--r', 'LineWidth',1.5, 'DisplayName','Delta 1st order');
            yline(L2_d2, ':k',  'LineWidth',1.5, 'DisplayName','Delta 2nd order');
            hold off;
            xlabel('GH quadrature order K');
            ylabel('L2 error vs MC Var (dB²)');
            title(sprintf('Var convergence — %s | %s | %s', plbl, tag, ref_lbl), ...
                  'Interpreter','none');
            legend('Location','best');
            grid on;
            saveFigPNG(fig_L2, fullfile(fig_dir, sprintf('delta_L2_orders_%s', pname)));
            drawnow; close(fig_L2);
        end

        %% Fig: EX difference — MC_EX vs TL_nom (3 params in one figure)
        if have_mc
            fig_ex = figure('Position',[50 50 1200 380]);
            EX_GH_K10 = {EX_GH_freq(:,:,GH_MAX_ORDER), ...
                         EX_GH_zS(:,:,GH_MAX_ORDER),   ...
                         EX_GH_svp(:,:,GH_MAX_ORDER)};
            ex_titles = {'E[TL] diff: GH-10 vs MC (freq)', ...
                         'E[TL] diff: GH-10 vs MC (zS)', ...
                         'E[TL] diff: GH-10 vs MC (svp)'};
            dex_max = 0;
            for pi = 1:3
                dex = EX_GH_K10{pi} - MC_EX_ref;
                dex_max = max(dex_max, max(abs(dex(:))));
            end
            dex_max = max(dex_max, 0.1);
            for pi = 1:3
                ax = subplot(1,3,pi);
                dex = EX_GH_K10{pi} - MC_EX_ref;
                pcolor(ax, r_km, z_m, dex);
                shading(ax,'interp'); set(ax,'YDir','reverse');
                colormap(ax, redblue(256)); colorbar(ax);
                clim(ax, [-dex_max, dex_max]);
                xlabel(ax,'Range (km)'); ylabel(ax,'Depth (m)');
                title(ax, ex_titles{pi}, 'FontSize',9,'Interpreter','none');
                overlayBathymetry(ax, bathy_m, sc.maxDepth_m);
            end
            sgtitle(sprintf('E[TL] diff (GH-10 − MC) | %s | %s (dB)', sc.name, dist.name), ...
                    'Interpreter','none');
            saveFigPNG(fig_ex, fullfile(fig_dir,'EX_diff_GH10_vs_MC'));
            drawnow; close(fig_ex);
        end

        %% Fig: Nominal TL vs |MC_EX - TL_nom| (per param, 3 panels)
        if have_mc
            fig_exn = figure('Position',[50 50 1200 380]);
            exn_max = max(abs(MC_EX_ref(:) - TL_nom(:)));
            exn_max = max(exn_max, 0.1);
            for pi = 1:3
                ax = subplot(1,3,pi);
                pcolor(ax, r_km, z_m, abs(MC_EX_ref - TL_nom));
                shading(ax,'interp'); set(ax,'YDir','reverse');
                colormap(ax, hot(256)); colorbar(ax);
                clim(ax, [0, exn_max]);
                xlabel(ax,'Range (km)'); ylabel(ax,'Depth (m)');
                title(ax, param_lbls{pi}, 'FontSize',9,'Interpreter','none');
                overlayBathymetry(ax, bathy_m, sc.maxDepth_m);
            end
            sgtitle(sprintf('|MC E[TL] − TL_{nom}| | %s | %s (dB)', sc.name, dist.name), ...
                    'Interpreter','none');
            saveFigPNG(fig_exn, fullfile(fig_dir,'EX_diff_MC_nominal'));
            drawnow; close(fig_exn);
        end

        fprintf('  Figures saved: %s\n', fig_dir);
    end
end

fprintf('\n=== run_Delta complete. ===\n');
end
