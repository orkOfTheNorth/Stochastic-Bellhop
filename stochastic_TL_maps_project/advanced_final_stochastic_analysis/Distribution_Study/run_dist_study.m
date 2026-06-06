%% run_dist_study.m  –  Delta Method Distribution & Error-Level Research
%
% PURPOSE
% ───────
% Systematically compare Delta method accuracy under:
%   (a) Different input distributions: Uniform vs Normal
%   (b) Different error levels: 1%, 5%, 10% of nominal value
%   (c) Single-parameter vs combined perturbation
%   (d) Bound comparison: Chebyshev, LN3, Empirical MC
%
% KEY INSIGHT
% ───────────
% The Delta Jacobian (∂TL/∂θ) is independent of the distribution —
% it is computed at the nominal point and reflects local sensitivity.
% Only the VARIANCE assigned to each parameter changes across configs.
%
% Therefore:
%   Var_delta = J² × σ²     changes ONLY with σ², not distribution shape
%   E_delta   = TL_nom      unchanged across all configs
%
% The comparison reveals where Delta breaks down (large errors → nonlinearity)
% and whether the bound method matters more than the distribution assumption.
%
% CONFIGS (5 total):
%   1. Uniform 1%  : U[nom±1%]   σ = (2%·nom)/√12
%   2. Uniform 5%  : U[nom±5%]   σ = (10%·nom)/√12
%   3. Uniform 10% : U[nom±10%]  σ = (20%·nom)/√12
%   4. Normal  5%  : N(nom, (5%·nom)²)
%   5. Normal  10% : N(nom, (10%·nom)²)
%
% SUBSETS ANALYZED (4 of 7):
%   freq only | zS only | svp only | zS+freq+svp (full)
%
% MC runs per config: 4 subsets × 50 = 200 runs → 5 configs = 1000 total runs
% Estimated time: ~3-5 s/run → 50-83 min

clear; close all; clc; warning('off');
try, cd(fileparts(mfilename('fullpath'))); catch; end

base_dir=fullfile('..');
addpath(genpath(fullfile(base_dir,'Functions')));
addpath(genpath(fullfile(base_dir,'../../early_stochastic_methods/code/Functions')));
addpath(base_dir);

if ~exist('results','dir'), mkdir('results'); end
if ~exist('figures','dir'), mkdir('figures'); end
tl_cache_dir='TL_cache'; if ~exist(tl_cache_dir,'dir'), mkdir(tl_cache_dir); end
cache_dir=fullfile(base_dir,'cache');

%% ── NOMINAL PARAMETERS (baseline scenario: const_35, 50km) ──────────────
freq0=10000; zS0=5; maxR=50000; FOM=100;
geo=[0.989*1500,1.63,0.07];
N=50;   % MC samples per config

%% ── LOAD OR COMPUTE JACOBIANS ────────────────────────────────────────────
jac_file=fullfile('..','Delta_Method','results','delta_jacobians.mat');
if ~isfile(jac_file)
    error(['Run Delta_Method/run_delta.m first to generate Jacobians.\n' ...
           'File expected at: %s'],jac_file);
end
J = load(jac_file,'J_freq','J_zS','J_svp','TL_nom','r_grid','z_grid','r_km','z_m');
J_freq=J.J_freq; J_zS=J.J_zS; J_svp=J.J_svp;
TL_nom=J.TL_nom; r_km=J.r_km; z_m=J.z_m; r_grid=J.r_grid; z_grid=J.z_grid;
[Nz,Nr]=size(TL_nom);
max_depth=max(z_grid);
fprintf('Jacobians loaded from %s\n',jac_file);

%% ── DISTRIBUTION CONFIGS ─────────────────────────────────────────────────
% Each config has: name, dist_type ('uniform'/'normal'), error_level (fraction)
% σ per parameter is computed from error_level and nominal value.
%
% For SVP noise: nominal scale = 25°C surface temp → 1% = 0.25°C
% (same as in run_delta.m and run_MC.m)
configs = struct( ...
    'name',   {'Uniform 1%','Uniform 5%','Uniform 10%','Normal 5%','Normal 10%'}, ...
    'dist',   {'uniform',   'uniform',   'uniform',    'normal',   'normal'}, ...
    'level',  {0.01,        0.05,        0.10,         0.05,       0.10} );

%% ── SUBSET SELECTION FOR THIS STUDY ─────────────────────────────────────
study_subsets = [1, 2, 3, 7];  % freq-only, zS-only, svp-only, full 3-param
study_names   = {'freq','zS','svp','zS_freq_svp'};
study_labels  = {'Freq only','z_S only','SVP only','Full (all 3)'};
study_active  = {[1 0 0],[0 1 0],[0 0 1],[1 1 1]};

%% ── MAIN LOOP ────────────────────────────────────────────────────────────
THRESHOLDS=[0.70 0.80 0.90 0.95];
all_results = struct();

for cfg_idx=1:numel(configs)
    cfg=configs(cfg_idx);
    fprintf('\n========== CONFIG %d/5: %s ==========\n',cfg_idx,cfg.name);

    % Parameter standard deviations for this config
    if strcmp(cfg.dist,'uniform')
        % U[nom*(1-L), nom*(1+L)] → σ² = (2L·nom)²/12
        sig_freq = freq0 * cfg.level * 2 / sqrt(12);
        sig_zS   = zS0   * cfg.level * 2 / sqrt(12);
        sig_svp  = 25    * cfg.level * 2 / sqrt(12);   % 25°C nominal
    else
        % N(nom, (L·nom)²) → σ = L·nom
        sig_freq = freq0 * cfg.level;
        sig_zS   = zS0   * cfg.level;
        sig_svp  = 25    * cfg.level;
    end
    var_freq_c = sig_freq^2;
    var_zS_c   = sig_zS^2;
    var_svp_c  = sig_svp^2;

    fprintf('  σ_freq=%.1fHz  σ_zS=%.4fm  σ_svp=%.4f°C\n',sig_freq,sig_zS,sig_svp);

    % Delta variance maps (per subset)
    V_freq_c = J_freq.^2 * var_freq_c;
    V_zS_c   = J_zS.^2   * var_zS_c;
    V_svp_c  = J_svp.^2  * var_svp_c;
    Vpool_c  = {V_freq_c, V_zS_c, V_svp_c};

    cfg_res = struct();
    for sub_idx=1:numel(study_subsets)
        fl=study_active{sub_idx}; sn=study_names{sub_idx}; lbl=study_labels{sub_idx};
        fprintf('\n  --- Subset: %s ---\n',lbl);

        % Delta variance and Chebyshev bound
        Vs = fl(1)*Vpool_c{1}+fl(2)*Vpool_c{2}+fl(3)*Vpool_c{3};
        d=TL_nom-FOM; cheb_lb=zeros(Nz,Nr); mask=d>0;
        cheb_lb(mask)=max(0,1-Vs(mask)./(Vs(mask)+d(mask).^2));

        % MC samples
        rng(42+cfg_idx*100+sub_idx);
        if strcmp(cfg.dist,'uniform')
            bnd_freq=[freq0-freq0*cfg.level, freq0+freq0*cfg.level];
            bnd_zS  =[zS0-zS0*cfg.level,     zS0+zS0*cfg.level];
            bnd_svp =[-25*cfg.level,          25*cfg.level];
            lhs=lhsdesign(N,3);
            samp_freq=bnd_freq(1)+lhs(:,1)*diff(bnd_freq);
            samp_zS  =bnd_zS(1)+lhs(:,2)*diff(bnd_zS);
            samp_svp =bnd_svp(1)+lhs(:,3)*diff(bnd_svp);
        else  % normal — LHS-stratified via inverse CDF transform
            lhs_n = lhsdesign(N,3);
            samp_freq = freq0 + sig_freq * norminv(lhs_n(:,1));
            samp_zS   = zS0   + sig_zS  * norminv(lhs_n(:,2));
            samp_svp  = sig_svp * norminv(lhs_n(:,3));
        end

        % TL cache key includes config
        safe_name = regexprep(sprintf('%s_%s',cfg.name,sn),'[^A-Za-z0-9_]','_');
        tl_file   = fullfile(tl_cache_dir, sprintf('TL_%s.mat', safe_name));
        use_cache = false;
        runs_todo = 1:N;
        TL_all    = zeros(Nz, Nr, N);

        if isfile(tl_file)
            tc = load(tl_file, 'N_saved','cfg_name');
            cfg_ok = isfield(tc,'cfg_name') && strcmp(tc.cfg_name, cfg.name);
            if cfg_ok
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
            end
        end

        if ~use_cache
            m_tl = matfile(tl_file, 'Writable', true);
            if numel(runs_todo) == N
                m_tl.TL_save   = zeros(Nz, Nr, N, 'single');
                m_tl.runs_done = false(1, N);
                m_tl.cfg_name  = cfg.name;
                m_tl.N_saved   = int32(0);
            end
            rd_vec  = m_tl.runs_done;
            n_done  = sum(rd_vec);
            lbl_tag = sprintf('%s|%s', cfg.name(1:min(7,end)), sn);
            t_start = tic;
            for iter = 1:numel(runs_todo)
                i = runs_todo(iter);
                fi=freq0; if fl(1),fi=samp_freq(i);end
                zi=zS0;   if fl(2),zi=samp_zS(i);end
                si=0;     if fl(3),si=samp_svp(i);end
                sim={fi,maxR,zi,0,0,"summer","const_35",geo,FOM};
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
                    progressBar(n_done, N, lbl_tag, t_start);
                end
            end
            fprintf('  [ALL SAVED] %s  (N=%d, %.0fMB)\n', tl_file, N, Nz*Nr*N*4/1e6);
        end

        MC_EX=mean(TL_all,3); MC_Var=var(TL_all,0,3); MC_PrFOM=mean(TL_all>FOM,3);
        d2=MC_EX-FOM; mc_cheb=zeros(Nz,Nr); mask2=d2>0;
        mc_cheb(mask2)=max(0,1-MC_Var(mask2)./(MC_Var(mask2)+d2(mask2).^2));

        fprintf('  LN3...\n'); LN3_prob=zeros(Nz,Nr);
        for zi=1:Nz,for ri=1:Nr,LN3_prob(zi,ri)=ln3prob(squeeze(TL_all(zi,ri,:)),FOM);end,end

        cfg_res.(sn)=struct( ...
            'Delta_Var',Vs,'Delta_Cheb',cheb_lb, ...
            'MC_EX',MC_EX,'MC_Var',MC_Var,'MC_PrFOM',MC_PrFOM, ...
            'MC_Cheb',mc_cheb,'LN3_prob',LN3_prob, ...
            'sig_freq',sig_freq,'sig_zS',sig_zS,'sig_svp',sig_svp);
    end  % subset loop

    % Save config results
    res_file=fullfile('results',sprintf('dist_%s.mat',regexprep(cfg.name,'[^A-Za-z0-9]','_')));
    save(res_file,'cfg_res','r_km','z_m','r_grid','z_grid','FOM','N');
    all_results.(regexprep(cfg.name,'[^A-Za-z0-9]','_'))=cfg_res;
    fprintf('\n  Config %d saved: %s\n',cfg_idx,res_file);
end

%% ── FIGURES ──────────────────────────────────────────────────────────────
fprintf('\nGenerating figures...\n');
dst_dir = fullfile('figures','dist_study');
if ~exist(dst_dir,'dir'), mkdir(dst_dir); end
cfg_names = fieldnames(all_results);

methods_fn  = {'Delta_Cheb','MC_Cheb','MC_PrFOM','LN3_prob'};
methods_lbl = {'Delta Cheb','MC Cheb','MC Empirical','LN3'};
N_m = numel(methods_fn);  N_sub = numel(study_names);

% ── Fig A: Var[TL] scaling — all configs, full subset ─────────────────
fig_A = figure('Visible','off','Position',[50 50 1600 600]);
for ci=1:numel(cfg_names)
    cr=all_results.(cfg_names{ci}); if ~isfield(cr,'zS_freq_svp'),continue;end
    ax=subplot(2,ceil(numel(cfg_names)/2),ci);
    pcolor(r_km,z_m,cr.zS_freq_svp.Delta_Var); shading interp; set(ax,'YDir','reverse');
    colormap(ax,hot); colorbar;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('Var | %s | sf=%.0fHz',configs(ci).name,cr.zS_freq_svp.sig_freq),'FontSize',8);
end
sgtitle('Delta Var[TL] — Full Subset, All Distributions');
saveFig(fig_A, fullfile(dst_dir,'var_scaling')); close(fig_A);

% ── Fig B: MC − Delta Var (nonlinearity) ───────────────────────────────
fig_B = figure('Visible','off','Position',[80 80 1600 600]);
for ci=1:numel(cfg_names)
    cr=all_results.(cfg_names{ci}); if ~isfield(cr,'zS_freq_svp'),continue;end
    dv=cr.zS_freq_svp.MC_Var-cr.zS_freq_svp.Delta_Var;
    ax=subplot(2,ceil(numel(cfg_names)/2),ci);
    pcolor(r_km,z_m,dv); shading interp; set(ax,'YDir','reverse');
    mx=max(abs(dv(:))); if mx==0,mx=1;end
    colormap(ax,redblue(256)); colorbar; clim([-mx mx]);
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('MC-Delta Var | %s',configs(ci).name),'FontSize',8);
end
sgtitle('MC - Delta Var[TL] — Nonlinearity Error');
saveFig(fig_B, fullfile(dst_dir,'var_nonlinearity')); close(fig_B);

% ── Per-config composite: 4 subsets × 4 methods ───────────────────────
for ci=1:numel(cfg_names)
    cr=all_results.(cfg_names{ci});
    fw2 = max(1600, 380*N_m);
    fig_c = figure('Visible','off','Position',[50 50 fw2 max(700, 200*N_sub)]);
    for sub_idx=1:N_sub
        sn_s=study_names{sub_idx}; if ~isfield(cr,sn_s),continue;end
        sr=cr.(sn_s);
        for mi=1:N_m
            panel=(sub_idx-1)*N_m+mi;
            ax=subplot(N_sub,N_m,panel);
            if strcmp(methods_fn{mi},'Delta_Var')
                pcolor(r_km,z_m,sr.(methods_fn{mi})); shading interp; set(ax,'YDir','reverse');
                colormap(ax,hot); colorbar;
            else
                pcolor(r_km,z_m,sr.(methods_fn{mi})); shading interp; set(ax,'YDir','reverse');
                colormap(ax,parula(256)); clim([0 1]); colorbar;
            end
            if sub_idx==1, title(methods_lbl{mi},'FontSize',9); end
            if mi==1, ylabel(sprintf('%s\nDepth (m)',study_labels{sub_idx}),'FontSize',8);
            else,     ylabel('Depth (m)','FontSize',7); end
            if sub_idx==N_sub, xlabel('Range (km)','FontSize',8); end
        end
    end
    sgtitle(sprintf('P(TL>%ddB) | %s | Rows=Subsets, Cols=Methods',FOM,configs(ci).name));
    fname=regexprep(configs(ci).name,'[^A-Za-z0-9]','_');
    saveFig(fig_c, fullfile(dst_dir, sprintf('shadow_%s',fname))); close(fig_c);
end

% ── LN3 across all distributions × all subsets ────────────────────────
fw3 = max(1600, 300*numel(cfg_names));
fig_ln3 = figure('Visible','off','Position',[50 50 fw3 max(600,180*N_sub)]);
for sub_idx=1:N_sub
    sn_s=study_names{sub_idx};
    for ci=1:numel(cfg_names)
        cr=all_results.(cfg_names{ci}); if ~isfield(cr,sn_s),continue;end
        panel=(sub_idx-1)*numel(cfg_names)+ci;
        ax=subplot(N_sub,numel(cfg_names),panel);
        pcolor(r_km,z_m,cr.(sn_s).LN3_prob); shading interp; set(ax,'YDir','reverse');
        colormap(ax,parula(256)); clim([0 1]); colorbar;
        if sub_idx==1, title(configs(ci).name,'FontSize',8); end
        if ci==1, ylabel(sprintf('%s\nDepth (m)',study_labels{sub_idx}),'FontSize',8);
        else,     ylabel('Depth (m)','FontSize',7); end
        if sub_idx==N_sub, xlabel('Range (km)','FontSize',8); end
    end
end
sgtitle(sprintf('LN3 P(TL>%ddB) — All Distributions x All Subsets',FOM));
saveFig(fig_ln3, fullfile(dst_dir,'ln3_all_distributions')); close(fig_ln3);

% ── Normal vs Uniform (matched σ) ─────────────────────────────────────
if isfield(all_results,'Uniform_5_') && isfield(all_results,'Normal_5_')
    U=all_results.Uniform_5_.zS_freq_svp; No=all_results.Normal_5_.zS_freq_svp;
    pairs={U.Delta_Cheb,No.Delta_Cheb,'Delta Cheb'; ...
           U.MC_PrFOM,  No.MC_PrFOM,  'MC Empirical'; ...
           U.LN3_prob,  No.LN3_prob,  'LN3'};
    fig_D=figure('Visible','off','Position',[120 120 1600 700]);
    for k=1:size(pairs,1)
        ax=subplot(2,3,k);
        pcolor(r_km,z_m,pairs{k,1}); shading interp; set(ax,'YDir','reverse');
        colormap(ax,parula(256)); clim([0 1]); colorbar;
        title(sprintf('Uniform 5%% | %s',pairs{k,3}),'FontSize',9);
        xlabel('Range (km)'); ylabel('Depth (m)');
        ax=subplot(2,3,3+k);
        pcolor(r_km,z_m,pairs{k,2}); shading interp; set(ax,'YDir','reverse');
        colormap(ax,parula(256)); clim([0 1]); colorbar;
        title(sprintf('Normal 5%% | %s',pairs{k,3}),'FontSize',9);
        xlabel('Range (km)'); ylabel('Depth (m)');
    end
    sgtitle(sprintf('Normal vs Uniform 5%% — P(TL>%ddB) | Full 3-Param',FOM));
    saveFig(fig_D, fullfile(dst_dir,'normal_vs_uniform_5pct')); close(fig_D);
end

% ── Single vs combined — Uniform 1% ───────────────────────────────────
if isfield(all_results,'Uniform_1_')
    cr1=all_results.Uniform_1_;
    fig_E=figure('Visible','off','Position',[140 140 1600 800]);
    for k=1:N_sub
        sn_s=study_names{k}; if ~isfield(cr1,sn_s),continue;end; sr=cr1.(sn_s);
        ax=subplot(2,N_sub,k);
        pcolor(r_km,z_m,sr.Delta_Var); shading interp; set(ax,'YDir','reverse');
        colormap(ax,hot); colorbar;
        title(sprintf('Var | %s',study_labels{k}),'FontSize',9);
        xlabel('Range (km)'); ylabel('Depth (m)');
        ax=subplot(2,N_sub,N_sub+k);
        shadowThresholdMaps(ax,r_km,z_m,sr.MC_PrFOM, ...
            sprintf('Emp P | %s',study_labels{k}),THRESHOLDS);
    end
    sgtitle('Single vs Combined | Uniform 1% | Top=Var, Bottom=Empirical P(shadow)');
    saveFig(fig_E, fullfile(dst_dir,'single_vs_combined')); close(fig_E);
end

fprintf('\n=== Distribution study complete. Figures in figures/dist_study/ ===\n');
