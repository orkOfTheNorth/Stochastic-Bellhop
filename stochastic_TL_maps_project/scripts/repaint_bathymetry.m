%% repaint_bathymetry.m
% Repaints bathymetry overlay onto existing Delta + MC PNGs for slope
% scenarios only (upslope, downslope). Loads saved .mat results, redraws
% all figures with seafloor line + fill, overwrites the PNGs, then
% deletes this script.

try, cd(fileparts(mfilename('fullpath'))); catch; end
addpath(genpath('Shared_Utils'));
addpath(genpath('Bellhop'));

cfg = loadConfig();
[subset_names, subset_labels, active] = subsetDefs();
FOM        = cfg.nominal.FOM_dB;
THRESHOLDS = [0.80, 0.90, 0.95];
SLOPE_SCENARIOS = {'upslope', 'downslope'};

for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);
    if ~ismember(sc.name, SLOPE_SCENARIOS), continue; end

    bathy_m = bathymetryMaker(sc.bathy_type, sc.maxR_m);

    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);
        fprintf('Repainting %s | %s ...\n', sc.name, dist.name);

        %% ── DELTA figures ────────────────────────────────────────────────
        delta_res = fullfile('Methods','Delta', sc.name, dist.name, 'results');
        delta_fig = fullfile('Methods','Delta', sc.name, dist.name, 'figures');

        if isfolder(delta_res)
            % Load shared grid from first available subset mat
            jac_file = fullfile(delta_res, 'delta_jacobians.mat');
            if isfile(jac_file)
                J = load(jac_file, 'r_km','z_m');
                r_km = J.r_km;  z_m = J.z_m;
            else
                continue;
            end

            dz = z_m(2)-z_m(1);  dr = (r_km(2)-r_km(1))*1000;
            SZ_M=1; SR_M=100;
            sz_px=SZ_M/dz; sr_px=SR_M/dr;
            MIN_FEAT_Z=20; MIN_FEAT_R=500;
            lpf_cut_z=min(0.90,2*dz/MIN_FEAT_Z);
            lpf_cut_r=min(0.90,2*dr/MIN_FEAT_R);

            for s = 1:7
                sn  = subset_names{s};
                lbl = subset_labels{s};
                sf  = fullfile(delta_res, sprintf('delta_%s.mat', sn));
                if ~isfile(sf), continue; end
                D = load(sf, 'Cheb_lb_s');
                C = D.Cheb_lb_s;

                % cheb_
                figA = figure('Visible','off','Position',[50 50 700 480]);
                shadowCategoryMap(gca, r_km, z_m, C, ...
                    sprintf('Delta Chebyshev | %s | %s | %s', sc.name, dist.name, lbl), THRESHOLDS);
                overlayBathymetry(gca, bathy_m, sc.maxDepth_m);
                saveFigPNG(figA, fullfile(delta_fig, sprintf('cheb_%s', sn)));
                close(figA);

                C_lpf  = lpf2d(C, lpf_cut_z, lpf_cut_r);
                C_gaus = gaussBlur2d(C, sz_px, sr_px);
                lpf_lbl  = sprintf('LPF (z>%.0fm, r>%.0fm)', MIN_FEAT_Z, MIN_FEAT_R);
                gaus_lbl = sprintf('Gaussian (sz=%.0fm, sr=%.0fm)', SZ_M, SR_M);

                % smooth_ 3-panel
                figB = figure('Visible','off','Position',[50 50 1600 480]);
                ax1=subplot(1,3,1); shadowCategoryMap(ax1,r_km,z_m,C,     'Non-smoothed',THRESHOLDS); overlayBathymetry(ax1,bathy_m,sc.maxDepth_m);
                ax2=subplot(1,3,2); shadowCategoryMap(ax2,r_km,z_m,C_lpf, lpf_lbl,      THRESHOLDS); overlayBathymetry(ax2,bathy_m,sc.maxDepth_m);
                ax3=subplot(1,3,3); shadowCategoryMap(ax3,r_km,z_m,C_gaus,gaus_lbl,     THRESHOLDS); overlayBathymetry(ax3,bathy_m,sc.maxDepth_m);
                sgtitle(sprintf('Delta Smoothing | %s | %s | %s | FOM=%ddB', sc.name, dist.name, lbl, FOM));
                saveFigPNG(figB, fullfile(delta_fig, sprintf('smooth_%s', sn)));
                close(figB);

                % smooth_lpf_
                figL = figure('Visible','off','Position',[50 50 700 480]);
                shadowCategoryMap(gca,r_km,z_m,C_lpf, sprintf('Delta %s | %s | %s | %s',lpf_lbl,sc.name,dist.name,lbl),THRESHOLDS);
                overlayBathymetry(gca, bathy_m, sc.maxDepth_m);
                saveFigPNG(figL, fullfile(delta_fig, sprintf('smooth_lpf_%s', sn)));
                close(figL);

                % smooth_gauss_
                figG = figure('Visible','off','Position',[50 50 700 480]);
                shadowCategoryMap(gca,r_km,z_m,C_gaus,sprintf('Delta %s | %s | %s | %s',gaus_lbl,sc.name,dist.name,lbl),THRESHOLDS);
                overlayBathymetry(gca, bathy_m, sc.maxDepth_m);
                saveFigPNG(figG, fullfile(delta_fig, sprintf('smooth_gauss_%s', sn)));
                close(figG);
            end
            fprintf('  Delta done.\n');
        end

        %% ── MC figures ───────────────────────────────────────────────────
        mc_res = fullfile('Methods','MC', sc.name, dist.name, 'results');
        mc_fig = fullfile('Methods','MC', sc.name, dist.name, 'figures');

        if ~isfolder(mc_res), continue; end

        % Load any subset to get grid
        test_f = fullfile(mc_res, 'MC_zS.mat');
        if ~isfile(test_f), fprintf('  MC not yet run, skipping.\n'); continue; end
        G = load(test_f, 'r_km','z_m','N');
        r_km = G.r_km;  z_m = G.z_m;  N = G.N;
        max_depth = sc.maxDepth_m;

        all_stats = struct();
        for s = 1:7
            sn = subset_names{s};
            mf = fullfile(mc_res, sprintf('MC_%s.mat', sn));
            if ~isfile(mf), continue; end
            M = load(mf, 'MC_EX','MC_Var','MC_PrFOM','Cheb_lb','LN3_prob');
            all_stats.(sn) = M;
        end

        for s = 1:7
            sn  = subset_names{s};
            lbl = subset_labels{s};
            if ~isfield(all_stats, sn), continue; end
            st = all_stats.(sn);

            % TL map
            fig = figure('Visible','off','Position',[50 50 780 500]);
            plotTLmap(gca, r_km, z_m, st.MC_EX, FOM);
            overlayBathymetry(gca, bathy_m, max_depth);
            title(sprintf('MC E[TL] | %s | %s | %s | N=%d', sc.name, dist.name, lbl, N), 'Interpreter','none');
            saveFigPNG(fig, fullfile(mc_fig, sprintf('TL_%s', sn)));
            close(fig);

            % Chebyshev
            figC = figure('Visible','off','Position',[50 50 700 480]);
            shadowCategoryMap(gca,r_km,z_m,st.Cheb_lb, sprintf('Cheb | %s | %s | %s',sc.name,dist.name,lbl),THRESHOLDS);
            overlayBathymetry(gca, bathy_m, max_depth);
            saveFigPNG(figC, fullfile(mc_fig, sprintf('cheb_%s', sn)));
            close(figC);

            % Empirical
            figE = figure('Visible','off','Position',[50 50 700 480]);
            shadowCategoryMap(gca,r_km,z_m,st.MC_PrFOM,sprintf('Empirical | %s | %s | %s',sc.name,dist.name,lbl),THRESHOLDS);
            overlayBathymetry(gca, bathy_m, max_depth);
            saveFigPNG(figE, fullfile(mc_fig, sprintf('empirical_%s', sn)));
            close(figE);
        end

        % Combined EX / Var panels
        names_done = fieldnames(all_stats);
        if ~isempty(names_done)
            fw = max(1600, 200*7);
            ex_vals = cellfun(@(n) all_stats.(n).MC_EX, names_done, 'UniformOutput',false);
            ex_clim = [max(50,min(cellfun(@(m)min(m(:)),ex_vals))), min(150,max(cellfun(@(m)max(m(:)),ex_vals)))];

            fig_ex = figure('Visible','off','Position',[50 50 fw 420]);
            for s = 1:7
                sn = subset_names{s};
                if ~isfield(all_stats,sn), continue; end
                ax = subplot(1,7,s);
                pcolor(r_km,z_m,all_stats.(sn).MC_EX); shading interp; set(ax,'YDir','reverse');
                colormap(ax,jet); colorbar(ax); clim(ex_clim);
                hold(ax,'on');
                TL_sm = movmean(movmean(all_stats.(sn).MC_EX,20,2),20,1);
                contour(r_km,z_m,TL_sm,[FOM FOM],'w-','LineWidth',1.0);
                hold(ax,'off');
                overlayBathymetry(ax, bathy_m, max_depth);
                xlabel(ax,'Range (km)'); if s==1, ylabel(ax,'Depth (m)'); end
                title(ax, subset_labels{s}, 'FontSize',8,'Interpreter','none');
            end
            sgtitle(sprintf('MC E[TL] — %s | %s | N=%d | FOM=%ddB', sc.name, dist.name, N, FOM));
            saveFigPNG(fig_ex, fullfile(mc_fig,'combined_EX'));
            close(fig_ex);

            v_vals = cellfun(@(n) all_stats.(n).MC_Var, names_done, 'UniformOutput',false);
            v_max  = max(cellfun(@(m) max(m(:)), v_vals));
            fig_var = figure('Visible','off','Position',[50 50 fw 420]);
            for s = 1:7
                sn = subset_names{s};
                if ~isfield(all_stats,sn), continue; end
                ax = subplot(1,7,s);
                pcolor(r_km,z_m,all_stats.(sn).MC_Var); shading interp; set(ax,'YDir','reverse');
                colormap(ax,hot); colorbar(ax); clim([0,v_max]);
                overlayBathymetry(ax, bathy_m, max_depth);
                xlabel(ax,'Range (km)'); if s==1, ylabel(ax,'Depth (m)'); end
                title(ax, subset_labels{s}, 'FontSize',8,'Interpreter','none');
            end
            sgtitle(sprintf('MC Var[TL] — %s | %s | N=%d', sc.name, dist.name, N));
            saveFigPNG(fig_var, fullfile(mc_fig,'combined_Var'));
            close(fig_var);
        end

        fprintf('  MC done.\n');
    end
end

fprintf('\nAll slope figures repainted.\n');

% Self-delete
this_file = [mfilename('fullpath') '.m'];
clear functions;
delete(this_file);
