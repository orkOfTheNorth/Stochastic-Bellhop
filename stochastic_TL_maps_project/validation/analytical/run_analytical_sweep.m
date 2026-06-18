%% run_analytical_sweep.m — Waveguide frequency sweep from config.json
%
% Reads sweep_frequencies_Hz from config.json and for each frequency:
%   1. Runs the analytical normal-mode model (waveguide_analytical)
%   2. Runs Bellhop with IDENTICAL parameters:
%        — constant SVP at c=1500 m/s
%        — rigid bottom (geo = [10, 10, 0]: high impedance, no attenuation)
%        — same D, zS, maxR as the analytical model
%   3. Saves a 3-panel comparison PNG to analytical/figures/

SCRIPT_DIR = fileparts(mfilename('fullpath'));       % validation/analytical/
ROOT       = fileparts(fileparts(SCRIPT_DIR));       % project root
try; cd(SCRIPT_DIR); catch; end

addpath(genpath(fullfile(ROOT, 'core')));
addpath(fullfile(ROOT, 'binaries'));

cfg = loadConfig();

fig_dir = fullfile(SCRIPT_DIR, 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% ── Fixed waveguide parameters ────────────────────────────────────────────
c    = 1500;                     % sound speed (m/s) — constant SVP
D    = 2500;                     % water depth (m) — deep_water scenario (const_2500)
zS   = cfg.nominal.zS_m;        % 5 m
FOM  = cfg.nominal.FOM_dB;

% Rigid bottom: very high impedance (c_ratio=10, rho=10, alpha=0)
% → reflection coefficient ≈ +1 (Neumann), matching analytical model
geo_rigid = [10, 10, 0];

maxR_m     = 50000;              % 50 km range
bathy_type = 'const_2500';
cache_dir  = fullfile(ROOT, 'Cache', 'deep_water', 'bellhop_raw');
if ~exist(cache_dir, 'dir'), mkdir(cache_dir); end

% Thorp absorption (dB/km)
thorp = @(f_kHz) 0.11*f_kHz.^2./(1+f_kHz.^2) + 44*f_kHz.^2./(4100+f_kHz.^2) ...
                 + 3e-5*f_kHz.^2 + 0.003;

freqs = cfg.analytical.sweep_frequencies_Hz;



for fi = 1:numel(freqs)
    freq  = freqs(fi);
    f_kHz = freq / 1000;
    alpha = thorp(f_kHz);    % dB/km — same formula for both models

    fprintf('\n=== f = %d Hz ===\n', freq);

    %% Analytical normal-mode TL
    r_ana = linspace(500, maxR_m, 500);   % avoid r=0
    z_ana = linspace(0, D, 200);
    TL_ana = waveguide_analytical(r_ana, z_ana, freq, c, D, zS, alpha);

    %% Bellhop TL — SAME params as analytical model
    % SVP: const_1500 (matches c=1500 m/s)
    % Bottom: geo_rigid (matches Neumann rigid bottom)
    svp_str = sprintf('const_%.0f', c);   % 'const_1500'
    sim = {freq, maxR_m, zS, 0, 0, svp_str, bathy_type, geo_rigid, FOM};
    [TL_bhp, r_bhp, z_bhp] = bellhopCached(sim, cache_dir);
    r_bhp_km = r_bhp / 1000;

    %% Shared adaptive colour scale from Bellhop (no forced calibration)
    tl_lo = max(40,  prctile(TL_bhp(:), 2));
    tl_hi = min(160, prctile(TL_bhp(:), 98));
    clim_shared = [tl_lo, tl_hi];

    %% Figure — 3 panels: Analytical | Bellhop | Range profile at zS
    fig = figure('Position',[50 50 1400 500]);

    subplot(1,3,1);
    pcolor(r_ana/1000, z_ana, TL_ana);
    shading interp; set(gca,'YDir','reverse');
    colormap(jet); colorbar; clim(clim_shared);
    hold on;
    contour(r_ana/1000, z_ana, movmean(movmean(TL_ana,20,2),20,1), ...
            [FOM FOM], 'w-', 'LineWidth', 1.5);
    hold off;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('Analytical (modes)  f=%dHz', freq));

    subplot(1,3,2);
    pcolor(r_bhp_km, z_bhp, TL_bhp);
    shading interp; set(gca,'YDir','reverse');
    colormap(jet); colorbar; clim(clim_shared);
    hold on;
    contour(r_bhp_km, z_bhp, movmean(movmean(TL_bhp,20,2),20,1), ...
            [FOM FOM], 'w-', 'LineWidth', 1.5);
    hold off;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('Bellhop (rigid bty)  f=%dHz', freq));

    subplot(1,3,3);
    [~, zi_src_ana] = min(abs(z_ana - zS));
    [~, zi_src_bhp] = min(abs(z_bhp - zS));
    plot(r_ana/1000, TL_ana(zi_src_ana,:), 'b-',  'LineWidth', 1.5, 'DisplayName','Analytical');
    hold on;
    plot(r_bhp_km,   TL_bhp(zi_src_bhp,:), 'r--', 'LineWidth', 1.5, 'DisplayName','Bellhop');
    yline(FOM, 'k--', 'LineWidth', 1, 'DisplayName', sprintf('FOM=%ddB',FOM));
    hold off;
    set(gca,'YDir','reverse'); grid on; legend('Location','best');
    xlabel('Range (km)'); ylabel('TL (dB)');
    title(sprintf('Range profile at z_S=%.0fm', zS));

    sgtitle(sprintf('Waveguide vs Bellhop | f=%dHz  D=%gm  c=%dm/s  \\alpha=%.2fdB/km  Rigid bottom', ...
                    freq, D, c, alpha));

    out_path = fullfile(fig_dir, sprintf('waveguide_%dHz', freq));
    saveFigPNG(fig, out_path);
    drawnow; close(fig);
    fprintf('  Saved: waveguide_%dHz.png\n', freq);
end


fprintf('\n=== Analytical sweep complete. Figures in analytical/figures/ ===\n');
