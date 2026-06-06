%% run_MC.m  –  Monte Carlo UQ | Scenario: Upslope (550m→50m, 10km)
% N=50, ~3-5 s/run → ~18-29 min total.

clearvars -except N; close all; clc; warning('off');
try, cd(fileparts(mfilename('fullpath'))); catch; end
base_dir=fullfile('..','..');
addpath(genpath(fullfile(base_dir,'Functions')));
addpath(genpath(fullfile(base_dir,'../../early_stochastic_methods/code/Functions')));
addpath(base_dir);
if ~exist('results','dir'),mkdir('results');end
if ~exist('figures','dir'),mkdir('figures');end
tl_cache_dir='TL_cache'; if ~exist(tl_cache_dir,'dir'),mkdir(tl_cache_dir);end
cache_dir=fullfile(base_dir,'cache');

SCN_LABEL='Upslope (550m→50m, 10km)'; bathy_type='slope_550_50'; maxR=10000;
freq0=10000; zS0=5; FOM=100; geo=[0.989*1500,1.63,0.07];
if ~exist('N','var'),N=50;end; force_rerun=false;
freq_bnd=[freq0*0.99,freq0*1.01]; zS_bnd=[zS0*0.99,zS0*1.01]; svp_bnd=[-0.25,0.25];
subset_names  = {'zS','freq','svp','zS_freq','zS_svp','freq_svp','zS_freq_svp'};
subset_labels = {'z_S only','Freq only','SVP only','z_S+Freq','z_S+SVP','Freq+SVP','z_S+Freq+SVP'};
active={[0 1 0],[1 0 0],[0 0 1],[1 1 0],[0 1 1],[1 0 1],[1 1 1]};
rng(42); lhs=lhsdesign(N,3);
samp_freq=freq_bnd(1)+lhs(:,1)*diff(freq_bnd);
samp_zS=zS_bnd(1)+lhs(:,2)*diff(zS_bnd);
samp_svp=svp_bnd(1)+lhs(:,3)*diff(svp_bnd);

set(0,'DefaultFigureVisible','off');
sim_nom={freq0,maxR,zS0,0,0,"summer",bathy_type,geo,FOM};
[TL_dummy,r_grid,z_grid]=bellhopCached(sim_nom,cache_dir); close all;
r_km=r_grid/1000; z_m=z_grid; [Nz,Nr]=size(TL_dummy); max_depth=max(z_grid);

%% ── SVP PROFILE INSPECTION ────────────────────────────────────────────────
fig_svp=figure('Name','MC SVP Profiles','Visible','off','Position',[50 50 520 640]);
svp_nom=makeSVPNoise(0,max_depth); hold on;
for ii=1:N
    svp_i=makeSVPNoise(samp_svp(ii),max_depth);
    plot(svp_i(:,2),svp_i(:,1),'Color',[0.65 0.65 0.65],'LineWidth',0.6);
end
plot(svp_nom(:,2),svp_nom(:,1),'k-','LineWidth',2.5); hold off;
set(gca,'YDir','reverse'); xlabel('Sound Speed (m/s)'); ylabel('Depth (m)');
title(sprintf('SVP Profiles  N=%d  svp bnd=[%.2f,%.2f] deg C',N,svp_bnd(1),svp_bnd(2)));
grid on; saveFig(fig_svp,fullfile('figures','svp_profiles')); close(fig_svp);

fprintf('=== MC | %s | N=%d ===\n',SCN_LABEL,N);
all_stats=struct();
for s=1:7
    fl=active{s}; sn=subset_names{s};
    fprintf('\n=== Subset %d/7: [%s] ===\n',s,subset_labels{s});
    tl_file   = fullfile(tl_cache_dir, sprintf('TL_%s.mat', sn));
    use_cache = false;
    runs_todo = 1:N;
    TL_all    = zeros(Nz, Nr, N);

    if isfile(tl_file) && ~force_rerun
        tc = load(tl_file, 'N_saved','freq_bnd','zS_bnd','svp_bnd');
        bounds_ok = isfield(tc,'svp_bnd') && isequal(tc.freq_bnd,freq_bnd) && ...
                    isequal(tc.zS_bnd,zS_bnd) && isequal(tc.svp_bnd,svp_bnd);
        if bounds_ok
            if tc.N_saved >= N
                tmp    = load(tl_file,'TL_save');
                TL_all = double(tmp.TL_save(:,:,1:N));
                use_cache = true;
                fprintf('  [CACHE HIT] %d/%d runs loaded\n', N, tc.N_saved);
            elseif tc.N_saved > 0
                rd = load(tl_file,'TL_save','runs_done');
                if isfield(rd,'runs_done') && numel(rd.runs_done) >= N
                    runs_todo = find(~rd.runs_done(1:N));
                    TL_all    = double(rd.TL_save);
                    fprintf('  [RESUMING] %d/%d done, %d remaining\n', ...
                            tc.N_saved, N, numel(runs_todo));
                end
            end
        else
            fprintf('  [CACHE MISS] bounds changed — fresh run\n');
        end
    end
    if ~use_cache
        m_tl = matfile(tl_file, 'Writable', true);
        if numel(runs_todo) == N
            m_tl.TL_save   = zeros(Nz, Nr, N, 'single');
            m_tl.runs_done = false(1, N);
            m_tl.freq_bnd  = freq_bnd;
            m_tl.zS_bnd    = zS_bnd;
            m_tl.svp_bnd   = svp_bnd;
            m_tl.N_saved   = int32(0);
        end
        rd_vec = m_tl.runs_done;
        n_done = sum(rd_vec);
        t_start = tic;
        for iter = 1:numel(runs_todo)
            i = runs_todo(iter);
            fi=freq0;if fl(1),fi=samp_freq(i);end
            zi=zS0;if fl(2),zi=samp_zS(i);end
            si=0;if fl(3),si=samp_svp(i);end
            sim={fi,maxR,zi,0,0,"summer",bathy_type,geo,FOM};
            if fl(3),sim{6}="custom";
                [TL,~,~]=bellhopCached(sim,cache_dir,'CustomSVP',makeSVPNoise(si,max_depth));
            else,[TL,~,~]=bellhopCached(sim,cache_dir);end
            TL_all(:,:,i)       = TL;
            m_tl.TL_save(:,:,i) = single(TL);
            rd_vec(i)           = true;
            m_tl.runs_done      = rd_vec;
            n_done              = n_done + 1;
            m_tl.N_saved        = int32(n_done);
            close all;
            if mod(n_done,10)==0 || n_done==N
                progressBar(n_done, N, sn, t_start);
            end
        end
        fprintf('  [ALL SAVED] %s  (N=%d, %.0fMB)\n', tl_file, N, Nz*Nr*N*4/1e6);
    end
    MC_EX=mean(TL_all,3); MC_Var=var(TL_all,0,3); MC_PrFOM=mean(TL_all>FOM,3);
    d=MC_EX-FOM; lb=zeros(Nz,Nr); mask=d>0;
    lb(mask)=max(0,1-MC_Var(mask)./(MC_Var(mask)+d(mask).^2)); Cheb_lb=lb;
    fprintf('  LN3...\n'); LN3_prob=zeros(Nz,Nr);
    for zi=1:Nz,for ri=1:Nr,LN3_prob(zi,ri)=ln3prob(squeeze(TL_all(zi,ri,:)),FOM);end,end
    fname=fullfile('results',sprintf('MC_%s.mat',sn));
    save(fname,'MC_EX','MC_Var','MC_PrFOM','Cheb_lb','LN3_prob','TL_all', ...
         'r_grid','z_grid','r_km','z_m','FOM','N','freq_bnd','zS_bnd','svp_bnd','freq0','zS0','-v7.3');
    fprintf('  Saved %s\n',fname);
    all_stats.(sn)=struct('MC_EX',MC_EX,'MC_Var',MC_Var,'MC_PrFOM',MC_PrFOM, ...
                          'Cheb_lb',Cheb_lb,'LN3_prob',LN3_prob);
end
%% ── COMPOSITE SUMMARY FIGURES ────────────────────────────────────────────────
if isfile(fullfile('results','delta_jacobians.mat'))
    raw_ref=load(fullfile('results','delta_jacobians.mat'),'TL_nom'); raw_ref=raw_ref.TL_nom;
else, raw_ref=all_stats.(subset_names{7}).MC_EX; end
TL_ref=movmean(movmean(raw_ref,20,2),20,1);
N_s=7; fw=max(1400,200*N_s);

fig1=figure('Visible','off','Position',[50 50 fw 380]);
for s=1:N_s; sn=subset_names{s}; if ~isfield(all_stats,sn),continue;end
    ax=subplot(1,N_s,s); pcolor(r_km,z_m,all_stats.(sn).MC_EX); shading interp; set(ax,'YDir','reverse');
    colormap(ax,jet); colorbar; clim([50 150]);
    hold on; contour(r_km,z_m,TL_ref,[FOM FOM],'w-','LineWidth',0.8); hold off;
    xlabel('Range (km)'); if s==1,ylabel('Depth (m)');end; title(subset_labels{s},'FontSize',8);
end; sgtitle(sprintf('MC E[TL]  N=%d  |  %s',N,SCN_LABEL));
saveFig(fig1,fullfile('figures','MC_all_subsets_mean_TL')); close(fig1);

fig2=figure('Visible','off','Position',[50 50 fw 380]);
for s=1:N_s; sn=subset_names{s}; if ~isfield(all_stats,sn),continue;end
    ax=subplot(1,N_s,s); pcolor(r_km,z_m,all_stats.(sn).MC_Var); shading interp; set(ax,'YDir','reverse');
    colormap(ax,hot); colorbar;
    xlabel('Range (km)'); if s==1,ylabel('Depth (m)');end; title(subset_labels{s},'FontSize',8);
end; sgtitle(sprintf('MC Var[TL]  N=%d  |  %s',N,SCN_LABEL));
saveFig(fig2,fullfile('figures','MC_all_subsets_variance')); close(fig2);

fig3=figure('Visible','off','Position',[50 50 fw 380]);
for s=1:N_s; sn=subset_names{s}; if ~isfield(all_stats,sn),continue;end
    ax=subplot(1,N_s,s);
    shadowThresholdMaps(ax,r_km,z_m,all_stats.(sn).Cheb_lb,subset_labels{s},[0.70 0.80 0.90 0.95]);
end; sgtitle(sprintf('MC Chebyshev P(shadow)  N=%d  |  %s',N,SCN_LABEL));
saveFig(fig3,fullfile('figures','MC_all_subsets_cheb')); close(fig3);

fig4=figure('Visible','off','Position',[50 50 fw 380]);
for s=1:N_s; sn=subset_names{s}; if ~isfield(all_stats,sn),continue;end
    ax=subplot(1,N_s,s);
    shadowThresholdMaps(ax,r_km,z_m,all_stats.(sn).MC_PrFOM,subset_labels{s},[0.70 0.80 0.90 0.95]);
end; sgtitle(sprintf('MC Empirical P(shadow)  N=%d  |  %s',N,SCN_LABEL));
saveFig(fig4,fullfile('figures','MC_all_subsets_empirical')); close(fig4);

fig5=figure('Visible','off','Position',[50 50 fw 380]);
for s=1:N_s; sn=subset_names{s}; if ~isfield(all_stats,sn),continue;end
    ax=subplot(1,N_s,s);
    shadowThresholdMaps(ax,r_km,z_m,all_stats.(sn).LN3_prob,subset_labels{s},[0.70 0.80 0.90 0.95]);
end; sgtitle(sprintf('MC LN3 P(shadow)  N=%d  |  %s',N,SCN_LABEL));
saveFig(fig5,fullfile('figures','MC_all_subsets_ln3')); close(fig5);

set(0,'DefaultFigureVisible','on');
fprintf('\n=== MC complete: %s ===\n',SCN_LABEL);
