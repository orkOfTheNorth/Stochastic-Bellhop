%% run_comparison.m  –  Delta vs Monte Carlo Comparison  (2026-05-26)
%
% For each of the 7 parameter subsets, loads Delta and MC results and produces:
%   1. E[TL] difference map   (MC_EX − Delta_EX)
%   2. Var[TL] difference map (MC_Var − Delta_Var)
%   3. Shadow-zone boundary comparison (Chebyshev vs MC empirical vs LN3)
%      showing which regions each method declares as "95% shadow"
%
% Requires run_delta.m and run_MC.m to have been executed first.

clear; close all; clc;
try, cd(fileparts(mfilename('fullpath'))); catch; end
if ~exist('results','dir'), mkdir('results'); end
if ~exist('figures','dir'), mkdir('figures'); end

delta_dir = '../Delta_Method/results';
mc_dir    = '../Monte_Carlo/results';

subset_names  = {'zS','freq','temp','zS_freq','zS_temp','freq_temp','zS_freq_temp'};
subset_labels = {'z_S only','Freq only','Temp only', ...
                 'z_S + Freq','z_S + Temp','Freq + Temp','z_S + Freq + Temp'};

FOM = 100;   freq0 = 10000;   zS0 = 5;   % nominal params (not in delta .mat)

%% ── ALL-SUBSET: EX AND VAR DIFFERENCE MAPS ──────────────────────────────────
fig_EX  = figure('Name','Comparison – E[TL] Difference (MC−Delta) All Subsets', ...
                 'Position',[30 30 1700 900]);
fig_Var = figure('Name','Comparison – Var[TL] Difference (MC−Delta) All Subsets', ...
                 'Position',[60 60 1700 900]);
fig_Shd = figure('Name','Comparison – Shadow Zone Boundary (All Subsets)', ...
                 'Position',[90 90 1700 900]);

for s = 1:7
    sn  = subset_names{s};
    lbl = subset_labels{s};

    df = fullfile(delta_dir, sprintf('delta_%s.mat', sn));
    mf = fullfile(mc_dir,    sprintf('MC_%s.mat',    sn));

    if ~isfile(df) || ~isfile(mf)
        fprintf('Skipping subset "%s" (files not found).\n', lbl);
        continue;
    end

    D = load(df);   % TL_expected, Var_TL, Cheb_lb, r_km, z_m
    M = load(mf);   % MC_EX, MC_Var, MC_PrFOM, Cheb_lb

    r_km = D.r_km;   z_m = D.z_m;

    diff_EX  = M.MC_EX  - D.TL_expected;
    diff_Var = M.MC_Var - D.Var_TL;

    % Shadow zone masks (95% threshold)
    shd_delta_cheb = D.Cheb_lb_s >= 0.95;   % Delta Chebyshev 95%
    shd_mc_cheb    = M.Cheb_lb >= 0.95;   % MC   Chebyshev 95%
    shd_mc_emp     = M.MC_PrFOM >= 0.95;  % MC empirical   95%

    % For full subset: also compute LN3 shadow mask (on saved TL_all)
    shd_ln3 = false(size(D.TL_expected));
    if strcmp(sn,'zS_freq_temp') && isfield(M,'TL_all')
        fprintf('  Computing LN3 shadow mask for full subset (slow)...\n');
        TLa = M.TL_all;
        [Nz,Nr,~] = size(TLa);
        prob_ln3 = zeros(Nz,Nr);
        for zi = 1:Nz
            for ri = 1:Nr
                samps = double(squeeze(TLa(zi,ri,:)));
                prob_ln3(zi,ri) = ln3prob(samps, FOM);
            end
        end
        shd_ln3 = prob_ln3 >= 0.95;
        save(fullfile('results','LN3_shadow_full.mat'), ...
             'prob_ln3','shd_ln3','r_km','z_m','FOM');
    end

    % Save comparison data
    save(fullfile('results',sprintf('compare_%s.mat',sn)), ...
         'diff_EX','diff_Var','shd_delta_cheb','shd_mc_cheb','shd_mc_emp', ...
         'r_km','z_m','FOM','lbl');

    % ── E[TL] difference subplot ─────────────────────────────────────────
    figure(fig_EX);
    subplot(2,4,s);
    pcolor(r_km, z_m, diff_EX); shading interp; set(gca,'YDir','reverse');
    mx = max(abs(diff_EX(:))); if mx==0, mx=1; end
    colormap(gca,redblue(256)); colorbar; clim([-mx mx]);
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('ΔE[TL]: MC−Delta | %s',lbl),'FontSize',8);

    % ── Var difference subplot ────────────────────────────────────────────
    figure(fig_Var);
    subplot(2,4,s);
    pcolor(r_km, z_m, diff_Var); shading interp; set(gca,'YDir','reverse');
    mx = max(abs(diff_Var(:))); if mx==0, mx=1; end
    colormap(gca,redblue(256)); colorbar; clim([-mx mx]);
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('ΔVar[TL]: MC−Delta | %s',lbl),'FontSize',8);

    % ── Shadow comparison subplot ─────────────────────────────────────────
    figure(fig_Shd);
    subplot(2,4,s);
    % Encode: 0=none, 1=Delta only, 2=MC_emp only, 3=both, 4=Delta+MC_cheb
    shd_map = double(shd_delta_cheb) + 2*double(shd_mc_emp);
    imagesc(r_km, z_m, shd_map); set(gca,'YDir','reverse');
    colormap(gca, [0.9 0.9 0.9; 0.2 0.5 0.9; 0.9 0.3 0.2; 0.2 0.7 0.2]);
    clim([0 3]); colorbar('Ticks',[0,1,2,3], ...
        'TickLabels',{'None','Delta 95%','MC 95%','Both 95%'});
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('Shadow Zones 95%% | %s',lbl),'FontSize',8);

    fprintf('  Subset %d/7 done: %s\n', s, lbl);
end

% Add empty subplots for slot 8
figure(fig_EX);  subplot(2,4,8); axis off;
sgtitle('MC − Delta: E[TL] Difference Maps for All 7 Subsets  (red=MC higher, blue=Delta higher)');
saveFig(fig_EX, fullfile('figures','compare_EX_diff_all'));

figure(fig_Var); subplot(2,4,8); axis off;
sgtitle('MC − Delta: Var[TL] Difference Maps for All 7 Subsets  (red=MC higher, blue=Delta higher)');
saveFig(fig_Var, fullfile('figures','compare_Var_diff_all'));

figure(fig_Shd); subplot(2,4,8); axis off;
sgtitle(sprintf(['Shadow-Zone Comparison (P(TL>%ddB)≥95%%) | All 7 Subsets\n' ...
    'Blue=Delta-Cheb only | Red=MC-Empirical only | Green=Both agree'],FOM));
saveFig(fig_Shd, fullfile('figures','compare_shadow_all'));

%% ── FULL-SUBSET: DETAILED 4-METHOD SHADOW MAP ───────────────────────────────
D = load(fullfile(delta_dir,'delta_zS_freq_temp.mat'));
M = load(fullfile(mc_dir,   'MC_zS_freq_temp.mat'));
r_km = D.r_km;  z_m = D.z_m;

if isfile(fullfile('results','LN3_shadow_full.mat'))
    L = load(fullfile('results','LN3_shadow_full.mat'));
    shd_ln3   = L.shd_ln3;
    prob_ln3  = L.prob_ln3;
else
    shd_ln3  = false(size(D.TL_expected));
    prob_ln3 = zeros(size(D.TL_expected));
end

shd_delta = D.Cheb_lb_s >= 0.95;
shd_emp   = M.MC_PrFOM >= 0.95;
shd_mcheb = M.Cheb_lb  >= 0.95;

fig_full = figure('Name','Comparison – Full 3-Param: 4 Methods Shadow','Position',[100 100 1600 800]);

methods = {shd_delta, shd_mcheb, shd_emp, shd_ln3};
mtitles = {'Delta: Chebyshev 95%','MC: Chebyshev 95%','MC: Empirical 95%','MC: LN3 95%'};
for k = 1:4
    subplot(2,2,k);
    imagesc(r_km, z_m, double(methods{k})); set(gca,'YDir','reverse');
    colormap(gca,[0.85 0.85 0.85; 0.2 0.6 0.2]);
    colorbar('Ticks',[0.25 0.75],'TickLabels',{'Not Shadow','Shadow 95%'});
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('%s | Full 3-Param | f=%dHz zS=%.1fm FOM=%ddB', ...
          mtitles{k}, freq0, zS0, FOM));
end
sgtitle(sprintf(['Shadow Zone (P(TL>%ddB)≥95%%) – 4 Methods Compared | Full 3-Parameter Perturbation\n' ...
    'Green=shadow (guaranteed), Grey=not shadow'], FOM));
saveFig(fig_full, fullfile('figures','compare_shadow_4methods_full'));

%% ── FULL-SUBSET: DIFFERENCE MAPS (EX and VAR with interactive cursor) ────────
diff_EX_full  = M.MC_EX  - D.TL_expected;
diff_Var_full = M.MC_Var - D.Var_TL;

fig_diff = figure('Name','Comparison – Full 3-Param EX/Var Difference (interactive)', ...
                  'Position',[150 100 1100 800]);

ax1 = subplot(2,1,1);
pcolor(r_km, z_m, diff_EX_full); shading interp; set(gca,'YDir','reverse');
mx = max(abs(diff_EX_full(:))); if mx==0, mx=1; end
colormap(ax1, redblue(256)); colorbar; clim([-mx mx]);
xlabel('Range (km)'); ylabel('Depth (m)');
title(sprintf(['MC − Delta: E[TL] Difference [dB] | Full 3-Param\n' ...
    'f=%dHz±1%%  zS=%.1fm±1%%  Temp±1°C  |  N=%d MC runs'], freq0, zS0, M.N));

ax2 = subplot(2,1,2);
pcolor(r_km, z_m, diff_Var_full); shading interp; set(gca,'YDir','reverse');
mx = max(abs(diff_Var_full(:))); if mx==0, mx=1; end
colormap(ax2, redblue(256)); colorbar; clim([-mx mx]);
xlabel('Range (km)'); ylabel('Depth (m)');
title(sprintf(['MC − Delta: Var[TL] Difference | Full 3-Param\n' ...
    'Positive = MC estimates higher variance (nonlinearity error in Delta)']));

% Interactive cursor
dcm = datacursormode(fig_diff);
dcm.Enable = 'on';
dcm.DisplayStyle = 'window';
set(dcm,'UpdateFcn', @(~,evt) diffTip(evt,ax1,ax2,r_km,z_m, ...
    M.MC_EX,D.TL_expected,M.MC_Var,D.Var_TL,diff_EX_full,diff_Var_full));

saveFig(fig_diff, fullfile('figures','compare_EX_Var_full_interactive'));

%% ── FULL-SUBSET: 5-METHOD SHADOW COMPARISON INCLUDING RF ─────────────────────
% The RF shadow map is on a subsampled grid (stride=5) vs the full Bellhop grid,
% so the RF panel has lower spatial resolution but covers the same domain.
rf_file = fullfile('..','Decision_Tree','results','RF_results.mat');
if isfile(rf_file)
    rf_res = load(rf_file, 'RF_mc_avg', 'r_km_map', 'z_m_map');
    RF_shd = rf_res.RF_mc_avg >= 0.95;
    r_rf   = rf_res.r_km_map;   % subsampled range axis (km)
    z_rf   = rf_res.z_m_map;    % subsampled depth axis (m)

    fig5 = figure('Name','5-Method Shadow (incl RF)','Position',[120 120 1950 420]);

    all5_shd = {shd_delta, shd_mcheb, shd_emp, shd_ln3, RF_shd};
    all5_r   = {r_km, r_km, r_km, r_km, r_rf};
    all5_z   = {z_m,  z_m,  z_m,  z_m,  z_rf};
    all5_lbl = {'Delta Chebyshev', 'MC Chebyshev', 'MC Empirical', ...
                'MC LN3', 'RF MC-averaged'};

    for k = 1:5
        ax = subplot(1,5,k);
        imagesc(all5_r{k}, all5_z{k}, double(all5_shd{k}));
        set(ax,'YDir','reverse');
        colormap(ax, [0.88 0.88 0.88; 0.15 0.62 0.28]);   % grey=clear, green=shadow
        clim([0 1]);
        xlabel('Range (km)','FontSize',7);
        if k==1, ylabel('Depth (m)','FontSize',8); end
        title(sprintf('%s\n(95%% threshold)', all5_lbl{k}),'FontSize',8);
        if k==5
            cb = colorbar;
            cb.Ticks = [0.25 0.75];
            cb.TickLabels = {'Clear','Shadow'};
            cb.FontSize = 7;
        end
    end
    sgtitle(sprintf('Shadow Zone P(TL>%ddB)≥95%%  |  Full 3-Param  |  5 Methods', FOM), ...
        'FontSize',11);
    saveFig(fig5, fullfile('figures','compare_shadow_5methods_RF'));
    close(fig5);
    fprintf('Saved compare_shadow_5methods_RF\n');
else
    fprintf('Note: RF results not found — run run_DT.m to include RF in shadow comparison.\n');
end

%% ── FULL-SUBSET: CONTINUOUS PROBABILITY MAPS — ALL METHODS ──────────────────
% Complements the binary 5-method figure by showing the FULL probability
% gradient (0–1) for each method.  Reveals:
%   • How confidently each method declares a zone as shadow (≥95% threshold)
%   • Where methods agree on probability magnitude (not just ≥95% boundary)
%   • Spatial smoothness differences: Delta/RF are smoother; MC Emp is noisier
%
% Methods shown (2 rows × 3 panels, last panel blank):
%   Row 1: Delta Chebyshev bound  |  MC Chebyshev bound  |  MC Empirical
%   Row 2: MC LN3 fit             |  RF MC-averaged       |
%
% The white dashed contour marks the 95% threshold on each map.
% ─────────────────────────────────────────────────────────────────────────────
if isfile(rf_file)
    % Build LN3 continuous map (load from saved file if available)
    if ~exist('prob_ln3','var') || all(prob_ln3(:)==0)
        if isfile(fullfile('results','LN3_shadow_full.mat'))
            L2 = load(fullfile('results','LN3_shadow_full.mat'),'prob_ln3');
            prob_ln3 = L2.prob_ln3;
        end
    end

    % Continuous probability for each method
    prob_delta = D.Cheb_lb_s;          % Delta Chebyshev lower bound
    prob_mcheb = M.Cheb_lb;            % MC Chebyshev bound
    prob_emp   = M.MC_PrFOM;           % MC empirical fraction
    prob_rf    = rf_res.RF_mc_avg;     % RF MC-averaged (subsampled grid)

    all6_prob  = {prob_delta, prob_mcheb, prob_emp, prob_ln3, prob_rf};
    all6_r     = {r_km,       r_km,       r_km,     r_km,     r_rf};
    all6_z     = {z_m,        z_m,        z_m,      z_m,      z_rf};
    all6_lbl   = {'Delta Chebyshev', 'MC Chebyshev', 'MC Empirical', ...
                  'MC LN3', 'RF MC-averaged'};
    all6_note  = {'(conservative bound)', '(Cheb from MC Var)', ...
                  '(direct sample count)', '(3-param lognormal)', '(Random Forest)'};

    fig6 = figure('Name','Continuous P(shadow) — All Methods','Position',[130 50 1800 700]);

    for k = 1:5
        ax = subplot(2, 3, k);
        pcolor(all6_r{k}, all6_z{k}, all6_prob{k}); shading interp;
        set(ax,'YDir','reverse');
        colormap(ax, parula(256)); clim([0 1]);
        cb = colorbar; cb.Label.String = 'P(TL>FOM)'; cb.FontSize = 7;
        hold on;
        % White dashed contour at 95% threshold
        try
            contour(all6_r{k}, all6_z{k}, all6_prob{k}, [0.95 0.95], ...
                    'w--', 'LineWidth', 1.5);
        catch
        end
        hold off;
        xlabel('Range (km)','FontSize',7);
        if mod(k-1,3)==0, ylabel('Depth (m)','FontSize',8); end
        title({sprintf('\\bf%s', all6_lbl{k}), all6_note{k}}, ...
              'FontSize',8,'Interpreter','tex');
    end

    subplot(2, 3, 6); axis off;
    text(0.5, 0.7, {'White dashed contour = 95% threshold', '', ...
         'Colour shows CONTINUOUS P(shadow)', ...
         'not just the binary above/below 95%.', '', ...
         'Dark blue (P≈0) = definitely clear.', ...
         'Yellow (P≈1) = guaranteed shadow.', '', ...
         'Method differences are most visible', ...
         'in the 0.5–0.9 probability gradient.'}, ...
         'Units','norm','HorizontalAlignment','center', ...
         'FontSize',8,'VerticalAlignment','middle');

    sgtitle({sprintf('Continuous P(TL>%ddB) — All 5 Methods | Full 3-Parameter Uncertainty', FOM), ...
             'Colour = actual shadow probability (0=clear → 1=shadow)  |  White dashed = 95% boundary'}, ...
            'FontSize', 10);
    saveFig(fig6, fullfile('figures','compare_shadow_prob_all_methods'));
    close(fig6);
    fprintf('Saved compare_shadow_prob_all_methods\n');
end

fprintf('\n=== Comparison complete. Results in results/  Figures in figures/ ===\n');

%% ── LOCAL FUNCTIONS ──────────────────────────────────────────────────────────
function saveFig(fig, base_path)
    print(fig, base_path, '-dpng', '-r300');
    try
        exportgraphics(fig, [base_path '.pdf'], 'ContentType','image', 'Resolution',300);
    catch
    end
end

function txt = diffTip(evt,ax1,ax2,r_km,z_m,MC_EX,Delta_EX,MC_Var,Delta_Var,dEX,dVar)
    pos = get(evt,'Position');
    ca  = get(evt,'Target');
    clicked_ax = ca.Parent;
    [~,ri] = min(abs(r_km - pos(1)));
    [~,zi] = min(abs(z_m  - pos(2)));
    ri = max(1,min(ri,length(r_km)));
    zi = max(1,min(zi,length(z_m)));

    if clicked_ax == ax1
        txt = {sprintf('Range: %.2f km',  r_km(ri)), ...
               sprintf('Depth: %.1f m',   z_m(zi)), ...
               '── E[TL] ──────────────', ...
               sprintf('MC:    %.3f dB',  MC_EX(zi,ri)), ...
               sprintf('Delta: %.3f dB',  Delta_EX(zi,ri)), ...
               sprintf('Diff:  %.3f dB',  dEX(zi,ri))};
    else
        txt = {sprintf('Range: %.2f km',  r_km(ri)), ...
               sprintf('Depth: %.1f m',   z_m(zi)), ...
               '── Var[TL] ─────────────', ...
               sprintf('MC:    %.5f',     MC_Var(zi,ri)), ...
               sprintf('Delta: %.5f',     Delta_Var(zi,ri)), ...
               sprintf('Diff:  %.5f',     dVar(zi,ri))};
    end
end

function prob = ln3prob(samps, FOM)
% 3-parameter lognormal P(TL > FOM)
    samps = double(samps(:));
    gam   = min(samps) * 0.95;
    sh    = samps - gam;
    sh(sh<=0) = 1e-6;
    params = lognfit(sh);
    if FOM > gam
        prob = 1 - logncdf(FOM - gam, params(1), params(2));
    else
        prob = 1.0;
    end
end

function cmap = redblue(n)
% Red-white-blue diverging colormap (n steps)
    if nargin<1, n=256; end
    r1 = [linspace(0.2,1,n/2); linspace(0.2,1,n/2); ones(1,n/2)]';
    r2 = [ones(1,n/2); linspace(1,0.2,n/2); linspace(1,0.2,n/2)]';
    cmap = [r1; r2];
end
