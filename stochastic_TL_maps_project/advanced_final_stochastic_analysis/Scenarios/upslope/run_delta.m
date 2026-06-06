%% run_delta.m  –  Delta Method UQ | Scenario: Upslope
% Bathymetry: slope 550m → 50m over 10 km (shoaling)

clear; close all; clc; warning('off');
try, cd(fileparts(mfilename('fullpath'))); catch; end

base_dir=fullfile('..','..');
addpath(genpath(fullfile(base_dir,'Functions')));
addpath(genpath(fullfile(base_dir,'../../early_stochastic_methods/code/Functions')));
addpath(base_dir);

if ~exist('results','dir'), mkdir('results'); end
if ~exist('figures','dir'), mkdir('figures'); end
cache_dir=fullfile(base_dir,'cache');

SCN_LABEL='Upslope (550m→50m, 10km)';
bathy_type='slope_550_50'; maxR=10000;
freq0=10000; zS0=5; FOM=100;
geo=[0.989*1500,1.63,0.07];
freq_bnd=[freq0*0.99,freq0*1.01]; zS_bnd=[zS0*0.99,zS0*1.01]; svp_bnd=[-0.25,0.25];
var_freq=diff(freq_bnd)^2/12; var_zS=diff(zS_bnd)^2/12; var_svp=diff(svp_bnd)^2/12;
h_freq=50; h_zS=0.025; h_svp=0.125;

fprintf('=== DELTA | %s ===\n',SCN_LABEL);
set(0,'DefaultFigureVisible','off');

fprintf('[1/7] Nominal ...\n');
sim_nom={freq0,maxR,zS0,0,0,"summer",bathy_type,geo,FOM};
[TL_nom,r_grid,z_grid]=bellhopCached(sim_nom,cache_dir);
r_km=r_grid/1000; z_m=z_grid; max_depth=max(z_grid);

fprintf('[2/7] freq+h ...\n'); sim=sim_nom; sim{1}=freq0+h_freq;
[TL_fUp,~,~]=bellhopCached(sim,cache_dir);
fprintf('[3/7] freq-h ...\n'); sim{1}=freq0-h_freq;
[TL_fDn,~,~]=bellhopCached(sim,cache_dir);
fprintf('[4/7] zS+h ...\n'); sim=sim_nom; sim{3}=zS0+h_zS;
[TL_zUp,~,~]=bellhopCached(sim,cache_dir);
fprintf('[5/7] zS-h ...\n'); sim{3}=zS0-h_zS;
[TL_zDn,~,~]=bellhopCached(sim,cache_dir);
fprintf('[6/7] svp+h ...\n'); sim=sim_nom; sim{6}="custom";
[TL_sUp,~,~]=bellhopCached(sim,cache_dir,'CustomSVP',makeSVPNoise(+h_svp,max_depth));
fprintf('[7/7] svp-h ...\n');
[TL_sDn,~,~]=bellhopCached(sim,cache_dir,'CustomSVP',makeSVPNoise(-h_svp,max_depth));

%% ── SVP PERTURBATION INSPECTION ─────────────────────────────────────────────
fig_svp = figure('Name','Delta SVP Perturbations','Visible','off','Position',[50 50 520 640]);
svp_nom  = makeSVPNoise(0,       max_depth);
svp_plus = makeSVPNoise(+h_svp,  max_depth);
svp_neg  = makeSVPNoise(-h_svp,  max_depth);
hold on;
plot(svp_neg(:,2),  svp_neg(:,1),  'b-', 'LineWidth',1.5, 'DisplayName',sprintf('-%.3f deg C',h_svp));
plot(svp_nom(:,2),  svp_nom(:,1),  'k-', 'LineWidth',2.5, 'DisplayName','Nominal (0 deg C)');
plot(svp_plus(:,2), svp_plus(:,1), 'r-', 'LineWidth',1.5, 'DisplayName',sprintf('+%.3f deg C',h_svp));
hold off;
set(gca,'YDir','reverse'); grid on; legend('Location','southeast');
xlabel('Sound Speed (m/s)'); ylabel('Depth (m)');
title(sprintf('SVP Perturbations for FD Jacobian  h_{svp}=%.3f deg C',h_svp));
saveFig(fig_svp, fullfile('figures','svp_perturbations'));
close(fig_svp);

set(0,'DefaultFigureVisible','on'); close all;

J_freq=computeJacobian(TL_fUp,TL_fDn,h_freq);
J_zS=computeJacobian(TL_zUp,TL_zDn,h_zS);
J_svp=computeJacobian(TL_sUp,TL_sDn,h_svp);
V_freq=J_freq.^2*var_freq; V_zS=J_zS.^2*var_zS; V_svp=J_svp.^2*var_svp;
save(fullfile('results','delta_jacobians.mat'), ...
     'TL_nom','J_freq','J_zS','J_svp','V_freq','V_zS','V_svp', ...
     'r_grid','z_grid','r_km','z_m','FOM','freq0','zS0', ...
     'var_freq','var_zS','var_svp','freq_bnd','zS_bnd','svp_bnd');

subset_names  = {'zS','freq','svp','zS_freq','zS_svp','freq_svp','zS_freq_svp'};
subset_labels = {'z_S only','Freq only','SVP only','z_S+Freq','z_S+SVP','Freq+SVP','z_S+Freq+SVP'};
active={[0 1 0],[1 0 0],[0 0 1],[1 1 0],[0 1 1],[1 0 1],[1 1 1]};
Vpool={V_freq,V_zS,V_svp}; Var_delta=cell(1,7); Cheb_lb=cell(1,7);
for s=1:7
    fl=active{s}; Vs=fl(1)*Vpool{1}+fl(2)*Vpool{2}+fl(3)*Vpool{3}; Var_delta{s}=Vs;
    d=TL_nom-FOM; lb=zeros(size(TL_nom)); mask=d>0;
    lb(mask)=max(0,1-Vs(mask)./(Vs(mask)+d(mask).^2)); Cheb_lb{s}=lb;
    TL_expected=TL_nom; Var_TL=Vs; Cheb_lb_s=lb; Cheb_shadow95=lb>=0.95;
    save(fullfile('results',sprintf('delta_%s.mat',subset_names{s})), ...
         'TL_expected','Var_TL','Cheb_lb_s','Cheb_shadow95','r_grid','z_grid','r_km','z_m','FOM');
end

TL_contour=movmean(movmean(TL_nom,20,2),20,1); thresholds=[0.70 0.80 0.90 0.95];
N_s = 7;  fw = max(1400, 200*N_s);

fig4 = figure('Name','Delta – Mean TL all subsets','Visible','off','Position',[50 50 fw 380]);
for s = 1:N_s
    ax = subplot(1,N_s,s);
    pcolor(r_km,z_m,TL_nom); shading interp; set(ax,'YDir','reverse');
    colormap(ax,jet); colorbar; clim([50 150]);
    hold on; contour(r_km,z_m,TL_contour,[FOM FOM],'w-','LineWidth',0.8); hold off;
    xlabel('Range (km)'); if s==1, ylabel('Depth (m)'); end
    title(subset_labels{s},'FontSize',8);
end
sgtitle(sprintf('Delta E[TL]  (= TL_{nom}) | %s',SCN_LABEL));
saveFig(fig4,fullfile('figures','delta_all_subsets_mean_TL'));  close(fig4);

fig5 = figure('Name','Delta – Var all subsets','Visible','off','Position',[50 50 fw 380]);
for s = 1:N_s
    ax = subplot(1,N_s,s);
    pcolor(r_km,z_m,Var_delta{s}); shading interp; set(ax,'YDir','reverse');
    colormap(ax,hot); colorbar;
    xlabel('Range (km)'); if s==1, ylabel('Depth (m)'); end
    title(subset_labels{s},'FontSize',8);
end
sgtitle(sprintf('Delta Var[TL] — All 7 Subsets | %s',SCN_LABEL));
saveFig(fig5,fullfile('figures','delta_all_subsets_variance'));  close(fig5);

fig6 = figure('Name','Delta – Cheb shadow all subsets','Visible','off','Position',[50 50 fw 380]);
for s = 1:N_s
    ax = subplot(1,N_s,s);
    shadowThresholdMaps(ax,r_km,z_m,Cheb_lb{s},subset_labels{s},thresholds);
end
sgtitle(sprintf('Delta Chebyshev P(shadow)  FOM=%ddB  Thresholds: 70/80/90/95%% | %s',FOM,SCN_LABEL));
saveFig(fig6,fullfile('figures','delta_all_subsets_cheb'));  close(fig6);
fprintf('\n=== Delta complete: %s ===\n',SCN_LABEL);
