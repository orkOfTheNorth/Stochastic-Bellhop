%% run_analytical_delta_validation.m
%
% Ground-truth validation of the 1st-order Delta method using the
% analytical normal-mode waveguide solution (deep water, zS perturbation only).
%
% METHOD:
%   TL(r,z,zS) is computed analytically via waveguide_analytical.m.
%   The Jacobian dTL/dzS is derived analytically by differentiating the
%   mode shapes: d/dzS[sin(k_zm*zS)] = k_zm*cos(k_zm*zS).
%
%   1st-order Delta variance:  Var_delta(r,z) = (dTL/dzS)^2 * sigma_zS^2
%
%   "True" MC variance:        Var_MC(r,z) from N samples of zS ~ dist
%
% Comparison for three distributions: Normal 1%, 5%, 10%
%
% Output:
%   analytical/figures/delta_validation_<dist>.png  — 4-panel figure per dist
%   analytical/figures/delta_validation_summary.png — L1 summary bar chart
%   analytical/results/validation_<dist>.mat

clear; close all; clc; warning('off');
SCRIPT_DIR = fileparts(mfilename('fullpath'));       % validation/analytical/
ROOT       = fileparts(fileparts(SCRIPT_DIR));       % project root
try; cd(SCRIPT_DIR); catch; end
addpath(genpath(fullfile(ROOT, 'core')));
addpath(fullfile(ROOT, 'binaries'));

%% ── Parameters ──────────────────────────────────────────────────────────────
freq  = 1000;           % Hz  — 1kHz tractable in deep water (~3370 modes vs 33700 at 10kHz)
c     = 1484;           % m/s (deep water sound speed)
D     = 2500;           % m   (deep water depth, matches config)
zS0   = 100;            % m   (deeper source — 5m is near surface, sin(k_zm*5)≈0 for most modes)
alpha = 0.07;           % dB/km Thorp absorption at 1 kHz
FOM   = 100;            % dB

Nr    = 200;
Nz    = 100;
r_vec = linspace(100, 50000, Nr);    % [1 x Nr] m
z_vec = linspace(1,   D-1,   Nz)';  % [Nz x 1] m  (avoid 0 and D)

N_MC  = 2000;           % large N for clean ground truth
rng(42);

% Distributions: Normal 1%, 5%, 10%
dists = struct('name', {'Normal_1pct','Normal_5pct','Normal_10pct'}, ...
               'sig',  {zS0*0.01/2,   zS0*0.05/2,  zS0*0.10/2});

fig_dir = fullfile(SCRIPT_DIR, 'figures');
res_dir = fullfile(SCRIPT_DIR, 'results');
for d = {fig_dir, res_dir}
    if ~exist(d{1},'dir'), mkdir(d{1}); end
end

%% ── Analytical Jacobian helper ───────────────────────────────────────────────
% dTL/dzS at nominal zS0:
%   p(r,z,zS) = sum_m  A_m * sin(k_zm*zS) * sin(k_zm*z) * H(r)
%   dp/dzS    = sum_m  A_m * k_zm*cos(k_zm*zS) * sin(k_zm*z) * H(r)
%   dTL/dzS   = -20/ln(10) * Re(dp/dzS) / |p|   (chain rule on -20log10|p|)
%
% Both p and dp/dzS computed together in one mode loop.

fprintf('=== Analytical Delta Validation ===\n');
m_max_theory = floor(2*pi*freq/c * D/pi);
fprintf('Deep water: D=%.0fm  freq=%.0fHz  zS0=%.1fm  (~%d propagating modes)\n', D, freq, zS0, m_max_theory);

k     = 2*pi*freq/c;
alpha_npm = alpha / (20*log10(exp(1))*1000);

P_nom  = zeros(Nz, Nr);   % pressure field at zS0
dP_dzS = zeros(Nz, Nr);   % dP/dzS at zS0

for m = 0:40000
    k_zm   = (m + 0.5) * pi / D;
    inside = k^2 - k_zm^2;
    if inside <= 0, break; end
    k_xm   = sqrt(inside) - 1i*alpha_npm;

    mode_zS  = sin(k_zm * zS0);
    dmode_dzS = k_zm * cos(k_zm * zS0);
    mode_z   = sin(k_zm .* z_vec);          % [Nz x 1]
    prop     = exp(-1i*(k_xm.*r_vec - pi/4)) ./ sqrt(real(k_xm).*r_vec);  % [1 x Nr]

    Nm   = sqrt(pi/2) / D;
    term = Nm .* (mode_z * prop);           % [Nz x Nr]
    P_nom  = P_nom  + mode_zS   .* term;
    dP_dzS = dP_dzS + dmode_dzS .* term;
end

% TL and analytical Jacobian
TL_nom  = -20*log10(max(abs(P_nom), 1e-20));   % [Nz x Nr]
J_anal  = -20/log(10) .* real(dP_dzS) ./ max(abs(P_nom).^2, 1e-30) .* real(conj(P_nom));
% Simpler: dTL/dzS = -20/ln(10) * Re(P* dP/dzS) / |P|^2
J_anal2 = -20/log(10) .* real(conj(P_nom) .* dP_dzS) ./ max(abs(P_nom).^2, 1e-30);
J_anal  = J_anal2;

fprintf('Nominal TL computed. Max|J_anal|=%.3f dB/m\n', max(abs(J_anal(:))));

%% ── Numerical FD Jacobian — computed with same mode loop (not waveguide_analytical) ──
h_zS = 0.025;   % m
P_p = zeros(Nz, Nr);  P_m = zeros(Nz, Nr);
for m = 0:40000
    k_zm  = (m + 0.5) * pi / D;
    inside = k^2 - k_zm^2;
    if inside <= 0, break; end
    k_xm  = sqrt(inside) - 1i*alpha_npm;
    mode_z = sin(k_zm .* z_vec);
    prop   = exp(-1i*(k_xm.*r_vec - pi/4)) ./ sqrt(real(k_xm).*r_vec);
    Nm     = sqrt(pi/2) / D;
    term   = Nm .* (mode_z * prop);
    P_p    = P_p + sin(k_zm*(zS0+h_zS)) .* term;
    P_m    = P_m + sin(k_zm*(zS0-h_zS)) .* term;
end
TL_p = -20*log10(max(abs(P_p), 1e-20));
TL_m = -20*log10(max(abs(P_m), 1e-20));
J_fd = (TL_p - TL_m) / (2*h_zS);

J_err_pct = mean(abs(J_anal(:) - J_fd(:))) / max(mean(abs(J_anal(:))), 1e-10) * 100;
fprintf('FD Jacobian error vs analytical: %.2f%%\n', J_err_pct);

%% ── Per-distribution comparison ─────────────────────────────────────────────
L1_delta = zeros(1, numel(dists));
L1_fd    = zeros(1, numel(dists));
r_km     = r_vec / 1000;

for di = 1:numel(dists)
    d    = dists(di);
    sig  = d.sig;

    % Delta variance (analytical Jacobian)
    Var_delta_anal = J_anal.^2 * sig^2;   % [Nz x Nr]
    % Delta variance (FD Jacobian)
    Var_delta_fd   = J_fd.^2   * sig^2;

    % MC ground truth — sample zS, compute TL analytically, get variance
    zS_samples = zS0 + sig * randn(N_MC, 1);
    zS_samples = max(zS_samples, 0.5);    % keep above surface

    % Precompute mode basis [Nz x Nr x M] — reuse across samples
    % For each mode m, term_m = Nm * sin(k_zm*z) * prop  [Nz x Nr]
    % P(zS) = sum_m sin(k_zm*zS) * term_m
    % Compute sin(k_zm*zS_i) for all samples at once: [M x N_MC]
    % then TL = -20log10|sum_m sin_m(zS_i) * term_m|  — vectorized over samples
    fprintf('  [%s] Building mode basis...\n', d.name);
    mode_terms = {};   kzm_vec = [];
    for m = 0:40000
        k_zm  = (m+0.5)*pi/D;
        if k^2 - k_zm^2 <= 0, break; end
        k_xm  = sqrt(k^2-k_zm^2) - 1i*alpha_npm;
        Nm    = sqrt(pi/2)/D;
        mode_terms{end+1} = Nm .* (sin(k_zm.*z_vec) * ...
            (exp(-1i*(k_xm.*r_vec-pi/4)) ./ sqrt(real(k_xm).*r_vec)));
        kzm_vec(end+1) = k_zm;
    end
    M = numel(kzm_vec);
    % Stack: [Nz*Nr x M]
    basis = reshape(cell2mat(cellfun(@(t) t(:), mode_terms, 'UniformOutput',false)),Nz*Nr,M);
    % sin(k_zm * zS_i): [M x N_MC]
    sin_zS = sin(kzm_vec(:) * zS_samples(:)');
    % P for all samples: [Nz*Nr x N_MC]
    P_all  = basis * sin_zS;
    TL_mc  = reshape(-20*log10(max(abs(P_all),1e-20)), Nz, Nr, N_MC);
    fprintf('  [%s] MC done (%d modes, %d samples).\n', d.name, M, N_MC);
    Var_MC = var(TL_mc, 0, 3);   % [Nz x Nr]
    EX_MC  = mean(TL_mc, 3);

    % Mask interference nulls — exclude pixels where nominal TL > 120 dB
    null_mask = TL_nom(:) <= 120;
    mean_VM = max(mean(Var_MC(null_mask)), 1e-10);
    L1_delta(di) = mean(abs(Var_delta_anal(null_mask) - Var_MC(null_mask))) / mean_VM;
    L1_fd(di)    = mean(abs(Var_delta_fd(null_mask)   - Var_MC(null_mask))) / mean_VM;
    fprintf('  [%s] %.0f%% pixels used (TL<=120dB mask)\n', d.name, 100*mean(null_mask));

    fprintf('  [%s] L1_analytical=%.4f  L1_FD=%.4f\n', d.name, L1_delta(di), L1_fd(di));

    % Save results
    save(fullfile(res_dir, sprintf('validation_%s.mat', d.name)), ...
         'Var_delta_anal','Var_delta_fd','Var_MC','EX_MC','TL_nom', ...
         'J_anal','J_fd','r_km','z_vec','sig', '-v7.3');

    %% Figure: 4 panels
    fig = figure('Position',[50 50 1400 380]);
    clr_lim = prctile(Var_MC(:), 99);

    ax1 = subplot(1,4,1);
    pcolor(ax1, r_km, z_vec, Var_MC); shading(ax1,'interp');
    set(ax1,'YDir','reverse'); colorbar(ax1); clim(ax1,[0 clr_lim]);
    title(ax1, sprintf('Var_{MC}  (N=%d)\n%s', N_MC, d.name),'FontSize',8,'Interpreter','none');
    xlabel(ax1,'Range (km)'); ylabel(ax1,'Depth (m)');

    ax2 = subplot(1,4,2);
    pcolor(ax2, r_km, z_vec, Var_delta_anal); shading(ax2,'interp');
    set(ax2,'YDir','reverse'); colorbar(ax2); clim(ax2,[0 clr_lim]);
    title(ax2, sprintf('Var_{\\Delta} analytical\nL1=%.3f', L1_delta(di)),'FontSize',8);
    xlabel(ax2,'Range (km)');

    ax3 = subplot(1,4,3);
    diff_map = Var_delta_anal - Var_MC;
    mx = prctile(abs(diff_map(:)), 99);
    pcolor(ax3, r_km, z_vec, diff_map); shading(ax3,'interp');
    set(ax3,'YDir','reverse'); colormap(ax3, redblue(256)); clim(ax3,[-mx mx]); colorbar(ax3);
    title(ax3, 'Var_{\Delta,anal} - Var_{MC}','FontSize',8);
    xlabel(ax3,'Range (km)');

    ax4 = subplot(1,4,4);
    ok = null_mask & Var_MC(:) > prctile(Var_MC(null_mask),10);
    scatter(ax4, Var_MC(ok), Var_delta_anal(ok), 3, 'filled', ...
            'MarkerFaceAlpha',0.3,'MarkerFaceColor',[0 0.45 0.74]);
    hold(ax4,'on');
    lv = max(prctile(Var_MC(ok),99), prctile(Var_delta_anal(ok),99));
    plot(ax4,[0 lv],[0 lv],'r-','LineWidth',1.5);
    R2 = corr(Var_MC(ok), Var_delta_anal(ok))^2;
    xlabel(ax4,'Var_{MC}'); ylabel(ax4,'Var_{\Delta,anal}');
    title(ax4, sprintf('Scatter\nR^2=%.3f', R2),'FontSize',8);
    grid(ax4,'on'); axis(ax4,'equal'); hold(ax4,'off');

    sgtitle(sprintf('Analytical Delta Validation | %s | \\sigma_{zS}=%.3fm', d.name, sig), ...
            'FontSize',10);

    out = fullfile(fig_dir, sprintf('delta_validation_%s', d.name));
    saveas(fig, [out '.png']); close(fig);
    fprintf('  Saved: %s.png\n', out);
end

%% ── Summary bar chart ────────────────────────────────────────────────────────
fig = figure('Position',[50 50 560 380]);
x   = 1:numel(dists);
b1  = bar(x-0.2, L1_delta*100, 0.35, 'FaceColor',[0 0.45 0.74]); hold on;
b2  = bar(x+0.2, L1_fd*100,    0.35, 'FaceColor',[0.85 0.33 0.1]);
yline(10, 'k--', '10%', 'LineWidth',1.2,'LabelHorizontalAlignment','left');
set(gca,'XTick',x,'XTickLabel',{dists.name},'XTickLabelRotation',15);
ylabel('Relative L1 error (%)');
title('1st-order Delta accuracy vs analytical ground truth','FontSize',9);
legend([b1 b2], 'Analytical Jacobian','FD Jacobian','Location','northwest');
grid on;

fprintf('\nSummary (Relative L1%%):\n');
fprintf('%-15s  Analytical  FD\n','Distribution');
for di=1:numel(dists)
    fprintf('%-15s  %7.2f%%  %7.2f%%\n', dists(di).name, L1_delta(di)*100, L1_fd(di)*100);
end

out = fullfile(fig_dir, 'delta_validation_summary');
saveas(fig, [out '.png']); close(fig);
fprintf('Saved: %s.png\n', out);
fprintf('\n=== run_analytical_delta_validation.m complete ===\n');
