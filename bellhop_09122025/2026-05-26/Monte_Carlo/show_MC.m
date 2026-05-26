%% show_MC.m  –  Interactive MC Viewer  (2026-05-26)
%
% Loads the FULL 3-parameter MC result and displays an interactive map.
%
% HOVER  → tooltip shows: E[TL], Var[TL], Chebyshev %, MC empirical %, LN3 %
% CLICK  → opens a separate figure with:
%             • 20-bar histogram of the 200 TL samples at that point
%             • 3-parameter Log-Normal (LN3) PDF curve overlaid
%             • Text panel with all statistics
%
% LN3 estimation (per grid point, on-demand):
%   threshold  γ  = min(samples) × 0.95   (shift below minimum)
%   then fit 2-param lognormal to (TL − γ)  →  μ_ln, σ_ln
%   P(TL > FOM) = 1 − logncdf(FOM − γ, μ_ln, σ_ln)

clear; close all; clc;
try, cd(fileparts(mfilename('fullpath'))); catch; end

%% ── LOAD DATA ────────────────────────────────────────────────────────────────
fname = fullfile('results','MC_zS_freq_temp.mat');
if ~isfile(fname)
    error('Run run_MC.m first to generate: %s', fname);
end
load(fname);  % loads MC_EX, MC_Var, MC_PrFOM, Cheb_lb, TL_all, r_km, z_m, FOM, N
[Nz, Nr, ~] = size(TL_all);

%% ── MAIN FIGURE: E[TL] MAP ───────────────────────────────────────────────────
fig_main = figure('Name', ...
    sprintf('MC Interactive Viewer | Full 3-Param | N=%d | Click for histogram',N), ...
    'Position',[50 50 900 650]);

ax = axes('Parent', fig_main);
pcolor(ax, r_km, z_m, MC_EX); shading(ax,'interp');
set(ax,'YDir','reverse');
colormap(ax, jet); cb = colorbar(ax); clim(ax,[50 150]);
ylabel(cb,'TL (dB)');
xlabel(ax,'Range (km)'); ylabel(ax,'Depth (m)');
title(ax, sprintf(['MC Expected TL | All 3 Params | N=%d\n' ...
    'freq±1%%  zS±1%%  Temp±1°C  |  f=%dHz, zS=%.1fm, FOM=%ddB\n' ...
    'Click any point for distribution details'], N, freq0, zS0, FOM), ...
    'FontSize',10);

% Panel showing last-clicked statistics
info_ax = axes('Parent', fig_main, 'Position',[0.75 0.02 0.23 0.30]);
axis(info_ax,'off');
info_txt = text(info_ax, 0.05, 0.95, 'Click a point...', ...
    'Units','normalized','VerticalAlignment','top','FontSize',9, ...
    'FontName','FixedWidth');

%% ── HOVER: show tooltip via data cursor ──────────────────────────────────────
dcm = datacursormode(fig_main);
dcm.Enable = 'on';
dcm.DisplayStyle = 'window';
set(dcm,'UpdateFcn', @(~,evt) hoverTip(evt, ax, r_km, z_m, ...
    MC_EX, MC_Var, MC_PrFOM, Cheb_lb, TL_all, FOM));

%% ── CLICK: open histogram figure ─────────────────────────────────────────────
set(fig_main,'WindowButtonDownFcn', @(~,~) onClickMap(ax, r_km, z_m, ...
    MC_EX, MC_Var, MC_PrFOM, Cheb_lb, TL_all, FOM, N, info_txt, freq0, zS0));

fprintf('Interactive viewer ready. Click on the map.\n');

%% ═══════════════════════════════════════════════════════════════════════════
%%  CALLBACK FUNCTIONS
%% ═══════════════════════════════════════════════════════════════════════════

function txt = hoverTip(evt, ax, r_km, z_m, MC_EX, MC_Var, MC_PrFOM, Cheb_lb, TL_all, FOM)
% Datacursor tooltip on hover/click
    pos = evt.Position;
    cp  = get(ax,'CurrentPoint');
    r_q = cp(1,1);  z_q = cp(1,2);
    [~,ri] = min(abs(r_km - r_q));
    [~,zi] = min(abs(z_m  - z_q));
    ri = max(1,min(ri,size(MC_EX,2)));
    zi = max(1,min(zi,size(MC_EX,1)));

    ex_v  = MC_EX(zi,ri);
    var_v = MC_Var(zi,ri);
    pr_mc = MC_PrFOM(zi,ri)*100;
    pr_ch = Cheb_lb(zi,ri)*100;

    % LN3 probability (on-demand)
    samps = squeeze(TL_all(zi,ri,:));
    pr_ln = ln3prob(samps, FOM) * 100;

    txt = {sprintf('Range:   %.2f km',   r_km(ri)), ...
           sprintf('Depth:   %.1f m',    z_m(zi)), ...
           '──────────────────', ...
           sprintf('E[TL]    = %.2f dB', ex_v), ...
           sprintf('Var[TL]  = %.3f',    var_v), ...
           sprintf('σ[TL]    = %.3f dB', sqrt(var_v)), ...
           '──────────────────', ...
           sprintf('P(TL>FOM) Chebyshev = %.1f %%', pr_ch), ...
           sprintf('P(TL>FOM) MC Empir. = %.1f %%', pr_mc), ...
           sprintf('P(TL>FOM) LN3 fit   = %.1f %%', pr_ln)};
end

function onClickMap(ax, r_km, z_m, MC_EX, MC_Var, MC_PrFOM, Cheb_lb, TL_all, FOM, N, info_txt, freq0, zS0)
% Left-click opens histogram + LN3 figure
    if ~strcmp(get(gcf,'SelectionType'),'normal'), return; end
    cp = get(ax,'CurrentPoint');
    r_q = cp(1,1);  z_q = cp(1,2);

    % Check click is inside axes
    xl = get(ax,'XLim'); yl = get(ax,'YLim');
    if r_q < xl(1)||r_q > xl(2)||z_q < yl(1)||z_q > yl(2), return; end

    [~,ri] = min(abs(r_km - r_q));
    [~,zi] = min(abs(z_m  - z_q));
    ri = max(1,min(ri,size(MC_EX,2)));
    zi = max(1,min(zi,size(MC_EX,1)));

    samps  = double(squeeze(TL_all(zi,ri,:)));
    ex_v   = MC_EX(zi,ri);
    var_v  = MC_Var(zi,ri);
    pr_mc  = MC_PrFOM(zi,ri)*100;
    pr_ch  = Cheb_lb(zi,ri)*100;
    [pr_ln, gam, mu_ln, sig_ln] = ln3prob(samps, FOM);
    pr_ln_pct = pr_ln * 100;

    % Update info panel
    info_str = sprintf(['R=%.2fkm  Z=%.1fm\n' ...
        'E[TL]  = %.2f dB\n' ...
        'Var    = %.4f\n' ...
        'σ      = %.4f dB\n' ...
        'Cheb  = %.1f %%\n' ...
        'MC emp = %.1f %%\n' ...
        'LN3    = %.1f %%'], ...
        r_km(ri), z_m(zi), ex_v, var_v, sqrt(var_v), pr_ch, pr_mc, pr_ln_pct);
    set(info_txt,'String',info_str);
    drawnow;

    % ── Histogram + LN3 figure ──────────────────────────────────────────────
    hfig = figure('Name', sprintf('TL Distribution @ R=%.2fkm Z=%.1fm', r_km(ri),z_m(zi)), ...
                  'Position',[980 100 600 480]);

    % 20-bar histogram (normalized to density)
    histogram(samps, 20, 'Normalization','pdf', ...
              'FaceColor',[0.3 0.6 0.9],'EdgeColor','w','FaceAlpha',0.7);
    hold on;

    % LN3 PDF curve
    x_lo = min(samps)*0.98;  x_hi = max(samps)*1.02;
    x_range = linspace(x_lo, x_hi, 300);
    y_ln3 = lognpdf(x_range - gam, mu_ln, sig_ln);
    plot(x_range, y_ln3, 'r-','LineWidth',2.5,'DisplayName','LN3 fit');

    % FOM line
    xline(FOM,'k--','LineWidth',2,'DisplayName',sprintf('FOM=%ddB',FOM));

    legend('MC samples (20 bins)','LN3 PDF fit',sprintf('FOM=%ddB',FOM), ...
           'Location','best');
    xlabel('TL (dB)'); ylabel('Probability Density');
    title(sprintf(['MC Sample Distribution @ R=%.2f km, Depth=%.1f m\n' ...
        'N=%d  |  f=%dHz, zS=%.1fm  |  P(TL>FOM): Cheb=%.1f%%, MC=%.1f%%, LN3=%.1f%%'], ...
        r_km(ri), z_m(zi), N, freq0, zS0, pr_ch, pr_mc, pr_ln_pct));

    % Statistics text box
    stats_str = sprintf(['E[TL]    = %.3f dB\n' ...
                         'Var[TL]  = %.4f\n' ...
                         'σ[TL]    = %.4f dB\n' ...
                         'LN3: γ=%.3f μ=%.3f σ=%.3f\n' ...
                         'Chebyshev  P(TL>%ddB) = %.2f%%\n' ...
                         'MC Empir.  P(TL>%ddB) = %.2f%%\n' ...
                         'LN3 fit    P(TL>%ddB) = %.2f%%'], ...
                        ex_v, var_v, sqrt(var_v), gam, mu_ln, sig_ln, ...
                        FOM, pr_ch, FOM, pr_mc, FOM, pr_ln_pct);
    annotation(hfig,'textbox',[0.58 0.55 0.39 0.38], ...
        'String',stats_str,'FitBoxToText','off','BackgroundColor','w', ...
        'FontSize',8,'FontName','FixedWidth','EdgeColor','k');

    hold off;
    drawnow;
end

%% ── HELPER: 3-parameter lognormal fit and P(TL > FOM) ───────────────────────
function [prob, gam, mu_ln, sig_ln] = ln3prob(samps, FOM)
% Estimate LN3 parameters and compute P(TL > FOM).
% gamma = threshold (shift parameter)  estimated as 95% of minimum sample.
% Remaining parameters fit via lognfit on shifted data.
    samps = double(samps(:));
    gam = min(samps) * 0.95;           % threshold slightly below min
    shifted = samps - gam;
    shifted(shifted <= 0) = 1e-6;      % guard against numerical zero
    params = lognfit(shifted);         % [mu_ln, sigma_ln]
    mu_ln  = params(1);
    sig_ln = params(2);
    % P(TL > FOM) = 1 - logncdf(FOM - gamma, mu_ln, sig_ln)
    if FOM > gam
        prob = 1 - logncdf(FOM - gam, mu_ln, sig_ln);
    else
        prob = 1.0;   % FOM is below threshold, always exceeded
    end
end
