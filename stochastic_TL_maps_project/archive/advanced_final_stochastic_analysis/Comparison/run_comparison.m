%% run_comparison.m  –  Delta vs MC Comparison  (2026-05-26 | Baseline)
%
% Loads Delta + MC results for all 7 subsets and produces:
%   1. E[TL] difference maps (MC − Delta)
%   2. Var[TL] difference maps (MC − Delta)
%   3. Shadow-probability comparison at 70%, 80%, 90%, 95%
%      for Chebyshev (Delta), Chebyshev (MC), Empirical, LN3
%   4. Full-subset 4-method continuous probability maps

clear; close all; clc;
try, cd(fileparts(mfilename('fullpath'))); catch; end

base_dir = fullfile('..');
addpath(genpath(fullfile(base_dir, 'Functions')));

if ~exist('results','dir'), mkdir('results'); end
if ~exist('figures','dir'), mkdir('figures'); end

delta_dir = '../Delta_Method/results';
mc_dir    = '../Monte_Carlo/results';

subset_names  = {'zS','freq','svp','zS_freq','zS_svp','freq_svp','zS_freq_svp'};
subset_labels = {'z_S only','Freq only','SVP only', ...
                 'z_S+Freq','z_S+SVP','Freq+SVP','z_S+Freq+SVP'};

FOM = 100;  freq0 = 10000;  zS0 = 5;
THRESHOLDS = [0.70 0.80 0.90 0.95];

%% ── ALL SUBSETS: EX AND VAR DIFFERENCE ──────────────────────────────────────
fig_EX  = figure('Name','Compare E[TL] Diff','Position',[30 30 1700 900]);
fig_Var = figure('Name','Compare Var Diff',  'Position',[60 60 1700 900]);

for s = 1:7
    sn  = subset_names{s};
    lbl = subset_labels{s};

    df = fullfile(delta_dir, sprintf('delta_%s.mat', sn));
    mf = fullfile(mc_dir,    sprintf('MC_%s.mat',    sn));
    if ~isfile(df) || ~isfile(mf)
        fprintf('Skipping %s (files not found).\n', lbl); continue;
    end

    D = load(df);
    M = load(mf);
    r_km = D.r_km;  z_m = D.z_m;

    diff_EX  = M.MC_EX  - D.TL_expected;
    diff_Var = M.MC_Var - D.Var_TL;

    figure(fig_EX); subplot(2,4,s);
    pcolor(r_km, z_m, diff_EX); shading interp; set(gca,'YDir','reverse');
    mx = max(abs(diff_EX(:))); if mx==0, mx=1; end
    colormap(gca,redblue(256)); colorbar; clim([-mx mx]);
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('ΔE[TL]: MC−Delta | %s',lbl),'FontSize',8);

    figure(fig_Var); subplot(2,4,s);
    pcolor(r_km, z_m, diff_Var); shading interp; set(gca,'YDir','reverse');
    mx = max(abs(diff_Var(:))); if mx==0, mx=1; end
    colormap(gca,redblue(256)); colorbar; clim([-mx mx]);
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('ΔVar: MC−Delta | %s',lbl),'FontSize',8);

    save(fullfile('results',sprintf('compare_%s.mat',sn)), ...
         'diff_EX','diff_Var','r_km','z_m','FOM','lbl');
    fprintf('  %d/7 done: %s\n', s, lbl);
end

figure(fig_EX);  subplot(2,4,8); axis off;
sgtitle('MC − Delta: E[TL] Difference  (red=MC higher, blue=Delta higher)');
saveFig(fig_EX, fullfile('figures','compare_EX_diff_all'));

figure(fig_Var); subplot(2,4,8); axis off;
sgtitle('MC − Delta: Var[TL] Difference');
saveFig(fig_Var, fullfile('figures','compare_Var_diff_all'));

%% ── ALL SUBSETS: 4-METHOD SHADOW AT EACH THRESHOLD ──────────────────────────
for thr_idx = 1:numel(THRESHOLDS)
    thr = THRESHOLDS(thr_idx);
    fig_s = figure('Name',sprintf('Shadow %.0f%%',thr*100), ...
                   'Position',[90+thr_idx*30, 90+thr_idx*30, 1700, 900], 'Visible','off');
    for s = 1:7
        sn  = subset_names{s};
        lbl = subset_labels{s};
        df = fullfile(delta_dir, sprintf('delta_%s.mat', sn));
        mf = fullfile(mc_dir,    sprintf('MC_%s.mat',    sn));
        if ~isfile(df) || ~isfile(mf), continue; end
        D = load(df);  M = load(mf);

        % 4-category map: 0=none, 1=Delta-Cheb only, 2=MC-Emp only, 3=both
        shd_delta = D.Cheb_lb_s >= thr;
        shd_emp   = M.MC_PrFOM  >= thr;
        shd_map   = double(shd_delta) + 2*double(shd_emp);

        subplot(2,4,s);
        imagesc(D.r_km, D.z_m, shd_map); set(gca,'YDir','reverse');
        colormap(gca, [0.9 0.9 0.9; 0.2 0.5 0.9; 0.9 0.3 0.2; 0.2 0.7 0.2]);
        clim([0 3]);
        colorbar('Ticks',[0,1,2,3], ...
            'TickLabels',{'None','Delta Cheb','MC Emp','Both'});
        xlabel('Range (km)'); ylabel('Depth (m)');
        title(sprintf('Shadow %.0f%% | %s',thr*100,lbl),'FontSize',8);
    end
    subplot(2,4,8); axis off;
    sgtitle(sprintf('Shadow Zone P(TL>%ddB)≥%.0f%% | Blue=Delta Cheb | Red=MC Emp | Green=Both', ...
            FOM, thr*100));
    saveFig(fig_s, fullfile('figures',sprintf('compare_shadow_%.0fpct_all',thr*100)));
    close(fig_s);
end

%% ── FULL SUBSET: 4-METHOD CONTINUOUS PROBABILITY MAPS ───────────────────────
df_full = fullfile(delta_dir, 'delta_zS_freq_svp.mat');
mf_full = fullfile(mc_dir,    'MC_zS_freq_svp.mat');

if isfile(df_full) && isfile(mf_full)
    D = load(df_full);
    M = load(mf_full);
    r_km = D.r_km;  z_m = D.z_m;

    prob_delta = D.Cheb_lb_s;
    prob_mcheb = M.Cheb_lb;
    prob_emp   = M.MC_PrFOM;
    prob_ln3   = M.LN3_prob;

    all_prob  = {prob_delta, prob_mcheb, prob_emp, prob_ln3};
    all_lbl   = {'Delta Chebyshev','MC Chebyshev','MC Empirical','MC LN3'};
    all_note  = {'(conservative bound)','(Cheb from MC Var)', ...
                 '(direct count)','(3-param lognormal)'};

    fig6 = figure('Name','Full 3-Param Continuous P(shadow)','Position',[100 50 1700 700]);
    for k = 1:4
        ax = subplot(2,2,k);
        shadowThresholdMaps(ax, r_km, z_m, all_prob{k}, ...
            sprintf('\\bf%s\n%s',all_lbl{k},all_note{k}), THRESHOLDS);
    end
    sgtitle(sprintf(['Continuous P(TL>%ddB) — 4 Methods | Full 3-Param\n' ...
        'Contours at 70%% (grey) / 80%% (orange) / 90%% (red) / 95%% (blue)'],FOM));
    saveFig(fig6, fullfile('figures','compare_shadow_4methods_continuous'));

    % Binary 4-panel (95% only for side-by-side clarity)
    fig7 = figure('Name','Full: Binary Shadow 95%','Position',[150 100 1600 400]);
    for k = 1:4
        subplot(1,4,k);
        imagesc(r_km, z_m, double(all_prob{k} >= 0.95)); set(gca,'YDir','reverse');
        colormap(gca,[0.85 0.85 0.85; 0.15 0.62 0.28]);
        clim([0 1]);
        colorbar('Ticks',[0.25 0.75],'TickLabels',{'Clear','Shadow'});
        xlabel('Range (km)'); if k==1, ylabel('Depth (m)'); end
        title(sprintf('%s',all_lbl{k}),'FontSize',9);
    end
    sgtitle(sprintf('Shadow Zone P(TL>%ddB) ≥ 95%% | Full 3-Param',FOM));
    saveFig(fig7, fullfile('figures','compare_shadow_4methods_95pct'));
    fprintf('Saved full-subset 4-method figures.\n');
end

fprintf('\n=== Comparison complete. ===\n');
