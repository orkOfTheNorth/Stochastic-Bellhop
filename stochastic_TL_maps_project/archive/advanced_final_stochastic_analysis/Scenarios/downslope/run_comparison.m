%% run_comparison.m – Delta vs MC Comparison | Scenario: Downslope
clear; close all; clc;
try, cd(fileparts(mfilename('fullpath'))); catch; end

base_dir = fullfile('..', '..');
addpath(genpath(fullfile(base_dir, 'Functions')));

if ~exist('results','dir'), mkdir('results'); end
if ~exist('figures','dir'), mkdir('figures'); end

delta_dir = 'results';  mc_dir = 'results';
SCN_LABEL = 'Downslope (50m→550m, 10km)';
FOM = 100;
THRESHOLDS = [0.70 0.80 0.90 0.95];

subset_names  = {'zS','freq','svp','zS_freq','zS_svp','freq_svp','zS_freq_svp'};
subset_labels = {'z_S only','Freq only','SVP only','z_S+Freq','z_S+SVP','Freq+SVP','z_S+Freq+SVP'};

fig_EX  = figure('Name','Compare EX','Position',[30 30 1700 900]);
fig_Var = figure('Name','Compare Var','Position',[60 60 1700 900]);

for s=1:7
    sn=subset_names{s}; lbl=subset_labels{s};
    df=fullfile(delta_dir,sprintf('delta_%s.mat',sn));
    mf=fullfile(mc_dir,sprintf('MC_%s.mat',sn));
    if ~isfile(df)||~isfile(mf), fprintf('Skip %s\n',lbl); continue; end
    D=load(df); M=load(mf);
    diff_EX=M.MC_EX-D.TL_expected; diff_Var=M.MC_Var-D.Var_TL;

    figure(fig_EX); subplot(2,4,s);
    pcolor(D.r_km,D.z_m,diff_EX); shading interp; set(gca,'YDir','reverse');
    mx=max(abs(diff_EX(:))); if mx==0,mx=1; end
    colormap(gca,redblue(256)); colorbar; clim([-mx mx]);
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('ΔE[TL] | %s',lbl),'FontSize',8);

    figure(fig_Var); subplot(2,4,s);
    pcolor(D.r_km,D.z_m,diff_Var); shading interp; set(gca,'YDir','reverse');
    mx=max(abs(diff_Var(:))); if mx==0,mx=1; end
    colormap(gca,redblue(256)); colorbar; clim([-mx mx]);
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('ΔVar | %s',lbl),'FontSize',8);
    fprintf('  %d/7 done: %s\n',s,lbl);
end
figure(fig_EX); subplot(2,4,8); axis off;
sgtitle(sprintf('MC−Delta E[TL] Diff | %s',SCN_LABEL));
saveFig(fig_EX,fullfile('figures','compare_EX_diff_all'));
figure(fig_Var); subplot(2,4,8); axis off;
sgtitle(sprintf('MC−Delta Var Diff | %s',SCN_LABEL));
saveFig(fig_Var,fullfile('figures','compare_Var_diff_all'));

% Multi-threshold shadow maps
for thr_idx=1:numel(THRESHOLDS)
    thr=THRESHOLDS(thr_idx);
    fig_s=figure('Visible','off','Position',[90 90 1700 900]);
    for s=1:7
        sn=subset_names{s}; lbl=subset_labels{s};
        df=fullfile(delta_dir,sprintf('delta_%s.mat',sn));
        mf=fullfile(mc_dir,sprintf('MC_%s.mat',sn));
        if ~isfile(df)||~isfile(mf), continue; end
        D=load(df); M=load(mf);
        shd_map=double(D.Cheb_lb_s>=thr)+2*double(M.MC_PrFOM>=thr);
        subplot(2,4,s);
        imagesc(D.r_km,D.z_m,shd_map); set(gca,'YDir','reverse');
        colormap(gca,[0.9 0.9 0.9; 0.2 0.5 0.9; 0.9 0.3 0.2; 0.2 0.7 0.2]);
        clim([0 3]); colorbar('Ticks',[0,1,2,3],'TickLabels',{'None','Delta','MC Emp','Both'});
        xlabel('Range (km)'); ylabel('Depth (m)');
        title(sprintf('Shadow %.0f%% | %s',thr*100,lbl),'FontSize',8);
    end
    subplot(2,4,8); axis off;
    sgtitle(sprintf('Shadow %.0f%% | %s',thr*100,SCN_LABEL));
    saveFig(fig_s,fullfile('figures',sprintf('compare_shadow_%.0fpct',thr*100)));
    close(fig_s);
end

% Full 3-param: 4-method continuous maps
mf_full=fullfile(mc_dir,'MC_zS_freq_svp.mat');
df_full=fullfile(delta_dir,'delta_zS_freq_svp.mat');
if isfile(mf_full)&&isfile(df_full)
    D=load(df_full); M=load(mf_full);
    all_prob={D.Cheb_lb_s, M.Cheb_lb, M.MC_PrFOM, M.LN3_prob};
    all_lbl={'Delta Cheb','MC Cheb','MC Empirical','MC LN3'};
    fig6=figure('Position',[100 50 1700 700]);
    for k=1:4
        ax=subplot(2,2,k);
        shadowThresholdMaps(ax,D.r_km,D.z_m,all_prob{k},all_lbl{k},THRESHOLDS);
    end
    sgtitle(sprintf('P(TL>%ddB) – 4 Methods | %s | 70/80/90/95%% contours',FOM,SCN_LABEL));
    saveFig(fig6,fullfile('figures','compare_shadow_4methods_continuous'));
end

fprintf('\n=== Comparison complete: %s ===\n',SCN_LABEL);
