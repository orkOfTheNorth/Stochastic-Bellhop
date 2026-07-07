%% generate_missing_figures.m
% Generates all figures that are currently missing from findings.tex:
%  1. KDE validation for shallow_water and upslope (Normal_5pct)
%  2. Variance accuracy summary for shallow_water, upslope, downslope (Normal_5pct)
%  3. LN3 MOM vs MLE scatter (deep_water, Normal_5pct, all params)
%  4. Bimodal pixel histogram (deep_water, Normal_5pct, SVP — CZ boundary pixel)
%
% Run from stochastic_TL_maps_project/ root:
%   cd stochastic_TL_maps_project
%   matlab -batch "addpath(genpath('.')); generate_missing_figures"

ROOT = fileparts(fileparts(mfilename('fullpath')));
try; cd(ROOT); catch; end
addpath(genpath(fullfile(ROOT,'core')));
addpath(fullfile(ROOT,'binaries'));
addpath(fullfile(ROOT,'pipeline'));
addpath(fullfile(ROOT,'scripts'));

cfg = loadConfig();
PARAMS = {'freq','zS','svp'};

fprintf('\n==============================\n');
fprintf(' generate_missing_figures.m\n');
fprintf('==============================\n\n');

%% ── 1. KDE validation for shallow_water and upslope ─────────────────────────
fprintf('--- 1. KDE validation (shallow_water, upslope) ---\n');
targets = {'shallow_water','upslope'};
for ti = 1:numel(targets)
    sc_name = targets{ti};
    dist_name = 'Normal_5pct';
    for pi = 1:numel(PARAMS)
        param = PARAMS{pi};
        out = fullfile('Methods','Comparison', sc_name, dist_name, 'figures', ...
                       sprintf('kde_validation_%s.png', param));
        if isfile(out)
            fprintf('  SKIP %s/%s/%s (exists)\n', sc_name, dist_name, param);
            continue;
        end
        fprintf('  Generating KDE validation: %s/%s/%s\n', sc_name, dist_name, param);
        % V bug fix: the raw TL cube is NOT in Methods/MC/.../results/MC_<param>.mat
        % (that file only holds aggregated stats like MC_EX/MC_Var) — it lives in
        % Cache/<sc>/<dist>/TL_<param>.mat under field 'TL_save'. This was why this
        % section always silently SKIPPED every combo with "no TL_cube field".
        mc_file  = fullfile('Methods','MC', sc_name, dist_name, 'results', ...
                           sprintf('MC_%s.mat', param));
        tl_cache = fullfile('Cache', sc_name, dist_name, sprintf('TL_%s.mat', param));
        if ~isfile(mc_file) || ~isfile(tl_cache)
            fprintf('    SKIP — MC_%s.mat or TL cache not found\n', param);
            continue;
        end
        mc = load(mc_file, 'r_km', 'z_m');
        FOM = getFOM(cfg, sc_name);

        % Compute KDE P(detect) per pixel — TL cube is [Nz Nr N]
        tl_data = load(tl_cache, 'TL_save');
        TL = double(tl_data.TL_save);

        % Ensure [Nz Nr N]
        if size(TL,3) < size(TL,1) && size(TL,3) < size(TL,2)
            % TL is [N Nz Nr] — permute
            TL = permute(TL,[2 3 1]);
        end
        [Nz, Nr, N] = size(TL);

        % MC empirical P(detect)
        P_MC = mean(TL <= FOM, 3);  % P(TL <= FOM) = detection probability

        % KDE P(detect): Gaussian kernel per pixel
        P_KDE = zeros(Nz, Nr);
        bw_vec = zeros(1, Nz*Nr);
        TL_flat = reshape(TL, Nz*Nr, N);
        for px = 1:Nz*Nr
            y = TL_flat(px,:)';
            sig = std(y);
            if sig < 1e-6; P_KDE(px) = double(mean(y) <= FOM); continue; end
            h = 1.06 * sig * N^(-0.2);   % Silverman bandwidth
            % P(TL <= FOM) = mean Phi((FOM - y_i)/h)
            P_KDE(px) = mean(normcdf((FOM - y) / h));
        end
        P_KDE = reshape(P_KDE, Nz, Nr);

        % Load range/depth axes
        if isfield(mc,'r_km'); r_km = mc.r_km; else; r_km = (0:Nr-1)/(Nr-1)*150; end
        if isfield(mc,'z_m');  z_m  = mc.z_m;  else; z_m  = (0:Nz-1)/(Nz-1)*3000; end

        % Bathymetry (range_m, depth_m) — used directly by overlayBathymetry,
        % which expects x-axis in km and y-axis in meters.
        % V bug fix: bathymetryMaker needs the scenario's bathy_type STRING
        % (e.g. 'const_50', 'slope_550_50'), not the scenario name itself.
        sc_idx  = find(strcmp({cfg.scenarios.name}, sc_name), 1);
        bathy_m = bathymetryMaker(cfg.scenarios(sc_idx).bathy_type, max(r_km)*1000);
        max_depth_m = max(z_m);

        % Plot
        fig = figure('Visible','off','Position',[10 10 1400 500]);
        clim_p = [0 1];

        % Panel 1: MC empirical
        ax1 = subplot(1,3,1);
        imagesc(ax1, r_km, z_m, P_MC, clim_p);
        set(ax1,'YDir','reverse'); colormap(ax1,'jet'); colorbar(ax1);
        xlabel(ax1,'Range (km)'); ylabel(ax1,'Depth (m)');
        title(ax1,sprintf('MC Empirical P(detect)\n%s | %s | %s', sc_name, dist_name, param),'Interpreter','none');
        overlayBathymetry(ax1, bathy_m, max_depth_m);

        % Panel 2: KDE
        ax2 = subplot(1,3,2);
        imagesc(ax2, r_km, z_m, P_KDE, clim_p);
        set(ax2,'YDir','reverse'); colormap(ax2,'jet'); colorbar(ax2);
        xlabel(ax2,'Range (km)');
        title(ax2,sprintf('KDE P(detect)\n(Silverman bw)'),'Interpreter','none');
        overlayBathymetry(ax2, bathy_m, max_depth_m);

        % Panel 3: Difference (symmetric-log color scale — see run_comparison.m
        % for the same pattern applied to variance-difference maps)
        ax3 = subplot(1,3,3);
        diff_map = P_KDE - P_MC;
        nz = abs(diff_map(diff_map ~= 0));
        if isempty(nz)
            lin_thresh = eps;
        else
            lin_thresh = max(prctile(nz, 20), eps);
        end
        Z = sign(diff_map) .* log10(1 + abs(diff_map)/lin_thresh);
        zmax = max(abs(Z(:)));
        if zmax == 0, zmax = 1; end
        imagesc(ax3, r_km, z_m, Z, [-zmax zmax]);
        set(ax3,'YDir','reverse');
        colormap(ax3, redblue(256));
        cb = colorbar(ax3);
        tick_z = linspace(-zmax, zmax, 7);
        tick_v = sign(tick_z) .* lin_thresh .* (10.^abs(tick_z) - 1);
        cb.Ticks = tick_z;
        cb.TickLabels = arrayfun(@(v) sprintf('%.2g', v), tick_v, 'UniformOutput', false);
        cb.Label.String = 'P(detect) diff (symlog)';
        xlabel(ax3,'Range (km)');
        title(ax3,sprintf('KDE - MC (abs err mean=%.3f)', mean(abs(diff_map(:)))),'Interpreter','none');
        overlayBathymetry(ax3, bathy_m, max_depth_m);

        sgtitle(sprintf('KDE Validation | %s | %s | %s', sc_name, dist_name, param), 'Interpreter','none');

        out_dir = fullfile('Methods','Comparison', sc_name, dist_name, 'figures');
        if ~exist(out_dir,'dir'); mkdir(out_dir); end
        % V bug fix: saveFigPNG appends '.png' itself (expects an extension-less
        % base path) — passing 'out' (which already ends in .png) produced
        % double-extensioned 'kde_validation_freq.png.png' files.
        [~, no_ext, ~] = fileparts(out);
        saveFigPNG(fig, fullfile(out_dir, no_ext));
        close(fig);
        fprintf('    Saved: %s\n', out);
    end
end

%% ── 2. Variance accuracy summary ─────────────────────────────────────────────
fprintf('\n--- 2. Variance accuracy summary ---\n');
compare_variance_spatial;   % runs all scenarios/distributions internally

%% ── 3. LN3 MOM vs MLE ───────────────────────────────────────────────────────
fprintf('\n--- 3. LN3 MOM vs MLE ---\n');
compare_ln3_mom_vs_mle;

%% ── 4. Bimodal pixel histogram ───────────────────────────────────────────────
fprintf('\n--- 4. Bimodal pixel histogram ---\n');
out_bim = fullfile('Methods','MC','deep_water','Normal_5pct','figures','bimodal_histogram_svp.png');
if isfile(out_bim)
    fprintf('  SKIP (exists)\n');
else
    % Load deep_water / Normal_5pct / SVP TL cube
    mc_file = fullfile('Methods','MC','deep_water','Normal_5pct','results','MC_svp.mat');
    if ~isfile(mc_file)
        fprintf('  SKIP — MC_svp.mat not found\n');
    else
        mc = load(mc_file);
        if isfield(mc,'TL_cube'); TL = mc.TL_cube;
        elseif isfield(mc,'TL_samples'); TL = mc.TL_samples;
        else; error('no TL field'); end
        if size(TL,3) < size(TL,1) && size(TL,3) < size(TL,2)
            TL = permute(TL,[2 3 1]);
        end
        [Nz, Nr, N] = size(TL);
        if isfield(mc,'r_km'); r_km = mc.r_km; else; r_km=(0:Nr-1)/(Nr-1)*150; end
        if isfield(mc,'z_m');  z_m  = mc.z_m;  else; z_m=(0:Nz-1)/(Nz-1)*2500; end

        % Find CZ boundary pixel: pixel with highest std among pixels with
        % mean TL in [80 110] dB (neither deep shadow nor near-source)
        Var_px = var(TL, 0, 3);   % [Nz Nr]
        Mean_px = mean(TL, 3);
        mask = (Mean_px >= 80) & (Mean_px <= 110);
        Var_masked = Var_px;
        Var_masked(~mask) = 0;
        [~, idx] = max(Var_masked(:));
        [iz, ir] = ind2sub([Nz Nr], idx);
        tl_vec = squeeze(TL(iz, ir, :));

        fig = figure('Visible','off','Position',[10 10 900 420]);

        % Left: histogram of 1000 samples at this pixel
        ax1 = subplot(1,2,1);
        histogram(ax1, tl_vec, 40, 'FaceColor',[0.2 0.5 0.8], 'EdgeColor','none', 'Normalization','pdf');
        hold(ax1,'on');
        % Overlay LN3 fit
        gamma = 0.95 * min(tl_vec);
        Y = tl_vec - gamma;
        Y(Y <= 0) = 1e-6;
        [mu_hat, sig_hat] = lognfit(Y);
        xx = linspace(min(tl_vec)-2, max(tl_vec)+2, 500);
        xx_shift = max(xx - gamma, 1e-6);
        pdf_ln3 = lognpdf(xx_shift, mu_hat, sig_hat);
        plot(ax1, xx, pdf_ln3, 'r-', 'LineWidth', 2);
        xlabel(ax1, 'TL (dB)'); ylabel(ax1, 'PDF');
        title(ax1, sprintf('CZ boundary pixel\nRange=%.1f km, Depth=%.0f m', r_km(ir), z_m(iz)), 'Interpreter','none');
        legend(ax1, 'MC samples (N=1000)', 'LN3 MLE fit', 'Location','northwest');
        grid(ax1,'on');

        % Right: TL map with marker
        mc_nom = fullfile('Methods','MC','deep_water','Normal_5pct','results','MC_svp.mat');
        ax2 = subplot(1,2,2);
        imagesc(ax2, r_km, z_m/1000, Mean_px, [60 140]);
        set(ax2,'YDir','reverse'); colormap(ax2,'jet'); colorbar(ax2);
        hold(ax2,'on');
        plot(ax2, r_km(ir), z_m(iz)/1000, 'w^','MarkerSize',12,'MarkerFaceColor','white','LineWidth',2);
        xlabel(ax2,'Range (km)'); ylabel(ax2,'Depth (km)');
        title(ax2,'Mean TL map — bimodal pixel marked','Interpreter','none');

        sgtitle(sprintf('Bimodal TL distribution at CZ boundary | SVP 5%% Normal\nStd=%.1f dB, showing LN3 fit failure', std(tl_vec)), 'Interpreter','none');

        out_dir = fileparts(out_bim);
        if ~exist(out_dir,'dir'); mkdir(out_dir); end
        saveFigPNG(fig, out_bim);
        close(fig);
        fprintf('  Saved: %s\n', out_bim);
    end
end

fprintf('\n=== generate_missing_figures.m DONE ===\n');
