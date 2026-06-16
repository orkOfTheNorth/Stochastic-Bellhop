%% save_frequency_maps.m  —  Analytical Acoustics: Multi-Frequency TL Maps
%
% Generates PNG figures at multiple frequencies for:
%   1. Method of Images — Shallow  (depth=10m, range=10m, zs=1m)
%   2. Method of Images — Deep Ocean  (depth=2500m, range=10km, zs=5m)
%   3. Waveguide Modal Solution  (D=2500m, range=10km, zs=5m)
%
% Frequency selection rationale:
%   Shallow  — lambda = [15 5 1.5 0.5 0.15] m  (depth=10m)
%   Deep     — lambda = [100 30 7.5 1.5 0.15] m (depth=2500m)
%   Waveguide— chosen to show 1, 2, 5, 10, 25 propagating modes
%
% Run this script from within analytical/ to regenerate all figures.

clear; close all; clc;
try, cd(fileparts(mfilename('fullpath'))); catch; end
if ~exist('figures','dir'), mkdir('figures'); end

c  = 1500;   % m/s (nominal sound speed)

%% ═══════════════════════════════════════════════════════════════════════════
%% 1.  METHOD OF IMAGES — SHALLOW
%% ═══════════════════════════════════════════════════════════════════════════
fprintf('=== 1/3  Method of Images: Shallow ===\n');

zs_sh      = 1;
max_r_sh   = 10;   % m
max_d_sh   = 10;   % m
freqs_sh   = [100, 300, 1000, 3000, 10000];
res        = 500;

fw = max(300, 220*numel(freqs_sh));
fig_sh = figure('Visible','off','Position',[50 50 fw 420]);
for fi = 1:numel(freqs_sh)
    f   = freqs_sh(fi);
    k   = 2*pi*f/c;
    lam = c/f;
    x_vec = linspace(0.01, max_r_sh, res);
    z_vec = linspace(0.01, max_d_sh, res);
    [X, Z] = meshgrid(x_vec, z_vec);
    R1 = sqrt(X.^2 + (Z - zs_sh).^2);
    R2 = sqrt(X.^2 + (Z + zs_sh).^2);
    P  = (exp(1i*k*R1)./R1) - (exp(1i*k*R2)./R2);
    TL = -10*log10(abs(P).^2 + eps);

    % Individual PNG
    fig_i = figure('Visible','off','Position',[50 50 560 440]);
    contourf(x_vec, z_vec, TL, 60, 'LineStyle','none');
    colormap(jet(60));  colorbar;  caxis([0 50]);
    set(gca,'YDir','reverse');  grid on;
    xlabel('Range (m)');  ylabel('Depth (m)');
    title(sprintf('Method of Images (Shallow)  f=%d Hz  \\lambda=%.2f m', f, lam));
    exportgraphics(fig_i, fullfile('figures', sprintf('MoI_shallow_f%dHz.png', f)), 'Resolution', 300);
    close(fig_i);
    fprintf('  f=%5d Hz  lambda=%.3f m  -> MoI_shallow_f%dHz.png\n', f, lam, f);

    % Composite panel
    ax = subplot(1, numel(freqs_sh), fi);
    contourf(ax, x_vec, z_vec, TL, 40, 'LineStyle','none');
    colormap(ax, jet(40));  caxis(ax, [0 50]);
    set(ax, 'YDir','reverse');
    xlabel('Range (m)');  if fi == 1, ylabel('Depth (m)'); end
    title(ax, sprintf('f=%d Hz\n\\lambda=%.2f m', f, lam), 'FontSize', 8);
end
sgtitle('Method of Images — Shallow  (zs=1m, depth=10m, 3D point source)');
exportgraphics(fig_sh, fullfile('figures', 'MoI_shallow_all_frequencies.png'), 'Resolution', 300);
close(fig_sh);
fprintf('  Composite -> MoI_shallow_all_frequencies.png\n\n');

%% ═══════════════════════════════════════════════════════════════════════════
%% 2.  METHOD OF IMAGES — DEEP OCEAN
%% ═══════════════════════════════════════════════════════════════════════════
fprintf('=== 2/3  Method of Images: Deep Ocean ===\n');

zs_d     = 5;
max_r_d  = 10000;   % m (10 km)
max_d_d  = 2500;    % m
freqs_d  = [15, 50, 200, 1000, 10000];
res_x = 600;  res_z = 400;

fw = max(300, 220*numel(freqs_d));
fig_d = figure('Visible','off','Position',[50 50 fw 420]);
for fi = 1:numel(freqs_d)
    f   = freqs_d(fi);
    k   = 2*pi*f/c;
    lam = c/f;
    x_vec = linspace(1, max_r_d, res_x);
    z_vec = linspace(0, max_d_d, res_z);
    [X, Z] = meshgrid(x_vec, z_vec);
    R1 = sqrt(X.^2 + (Z - zs_d).^2);
    R2 = sqrt(X.^2 + (Z + zs_d).^2);
    P  = (exp(1i*k*R1)./(R1+eps)) - (exp(1i*k*R2)./(R2+eps));
    TL = -20*log10(abs(P) + eps);

    % Individual PNG
    fig_i = figure('Visible','off','Position',[50 50 700 500]);
    imagesc(x_vec/1000, z_vec, TL);
    colormap(jet);  colorbar;  caxis([60 110]);
    set(gca,'YDir','reverse');
    xlabel('Range (km)');  ylabel('Depth (m)');
    title(sprintf('Method of Images (Deep Ocean)  f=%d Hz  \\lambda=%.1f m', f, lam));
    exportgraphics(fig_i, fullfile('figures', sprintf('MoI_deep_f%dHz.png', f)), 'Resolution', 300);
    close(fig_i);
    fprintf('  f=%5d Hz  lambda=%6.2f m  -> MoI_deep_f%dHz.png\n', f, lam, f);

    % Composite panel
    ax = subplot(1, numel(freqs_d), fi);
    imagesc(ax, x_vec/1000, z_vec, TL);
    colormap(ax, jet);  caxis(ax, [60 110]);
    set(ax, 'YDir','reverse');
    xlabel('Range (km)');  if fi == 1, ylabel('Depth (m)'); end
    title(ax, sprintf('f=%d Hz\n\\lambda=%.1f m', f, lam), 'FontSize', 8);
end
sgtitle('Method of Images — Deep Ocean  (zs=5m, D=2500m, 10km range)');
exportgraphics(fig_d, fullfile('figures', 'MoI_deep_all_frequencies.png'), 'Resolution', 300);
close(fig_d);
fprintf('  Composite -> MoI_deep_all_frequencies.png\n\n');

%% ═══════════════════════════════════════════════════════════════════════════
%% 3.  WAVEGUIDE MODAL SOLUTION
%% ═══════════════════════════════════════════════════════════════════════════
fprintf('=== 3/3  Waveguide Modal Solution ===\n');

zs_wg    = 5;
D        = 2500;
max_r_wg = 10000;
freqs_wg = [15, 30, 75, 150, 300];
res_x_wg = 500;  res_z_wg = 300;

fw = max(300, 220*numel(freqs_wg));
fig_wg = figure('Visible','off','Position',[50 50 fw 420]);
for fi = 1:numel(freqs_wg)
    f   = freqs_wg(fi);
    k   = 2*pi*f/c;
    x_vec = linspace(0.1, max_r_wg, res_x_wg);
    z_vec = linspace(0, D, res_z_wg);
    [X, Z] = meshgrid(x_vec, z_vec);
    P = zeros(size(X));
    n_prop = 0;
    M_limit = min(2000, ceil(2*D*f/c) + 30);   % adaptive: cover all propagating modes
    for m = 0:M_limit
        k_zm   = (m + 0.5)*pi/D;
        inside = k^2 - k_zm^2;
        if inside > 0
            k_xm   = sqrt(inside);
            n_prop = n_prop + 1;
        else
            k_xm = -1i*sqrt(abs(inside));
        end
        term = (1/(1i*k_xm*D)) .* sin(k_zm*zs_wg) .* sin(k_zm*Z) .* exp(-1i*k_xm*X);
        P = P + term;
    end
    TL = -20*log10(abs(P));

    % Individual PNG
    fig_i = figure('Visible','off','Position',[50 50 700 500]);
    imagesc(x_vec, z_vec, TL);
    colormap(jet);  colorbar;  caxis([20 100]);
    set(gca,'YDir','reverse');
    xlabel('Range (m)');  ylabel('Depth (m)');
    title(sprintf('Waveguide Modal  f=%d Hz  D=%dm  %d propagating modes', f, D, n_prop));
    exportgraphics(fig_i, fullfile('figures', sprintf('waveguide_modal_f%dHz.png', f)), 'Resolution', 300);
    close(fig_i);
    fprintf('  f=%3d Hz  M_limit=%4d  %3d prop. modes  -> waveguide_modal_f%dHz.png\n', ...
            f, M_limit, n_prop, f);

    % Composite panel
    ax = subplot(1, numel(freqs_wg), fi);
    imagesc(ax, x_vec, z_vec, TL);
    colormap(ax, jet);  caxis(ax, [20 100]);
    set(ax, 'YDir','reverse');
    xlabel('Range (m)');  if fi == 1, ylabel('Depth (m)'); end
    title(ax, sprintf('f=%d Hz\n%d modes', f, n_prop), 'FontSize', 8);
end
sgtitle(sprintf('Waveguide Modal Solution  (zs=%dm, D=%dm, range=%dkm)', ...
        zs_wg, D, max_r_wg/1000));
exportgraphics(fig_wg, fullfile('figures', 'waveguide_modal_all_frequencies.png'), 'Resolution', 300);
close(fig_wg);
fprintf('  Composite -> waveguide_modal_all_frequencies.png\n\n');

fprintf('=== All maps saved to analytical/figures/ ===\n');
fprintf('Individual PNGs: %d files\n', 2*numel(freqs_sh)+numel(freqs_d)+1+numel(freqs_wg)+1+1);
