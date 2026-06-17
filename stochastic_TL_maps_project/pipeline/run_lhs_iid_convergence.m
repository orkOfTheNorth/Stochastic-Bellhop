function run_lhs_iid_convergence(sc_target, dist_target)
%% run_lhs_iid_convergence  LHS+MoM vs IID+MLE convergence study.
%
% Runs N=100 Bellhop calls per strategy, sub-sampled at N_vec=[5 10 20 30 50 75 100].
% LHS + moment-matching LN3  (ln3Moments) -- stratified sampling
% IID + profile-MLE LN3      (ln3MLE)     -- i.i.d. sampling, MLE fit
%
% Both converge to the same answer at large N.  RMSE between them shrinks
% to zero, confirming both methods are unbiased.
% NOTE: final production P(detect) results come from IID+MLE (run_mc).
%
% Outputs: Methods/LHS_IID_Conv/<scen>/<dist>/conv_<param>.png

if nargin < 1, sc_target  = ''; end
if nargin < 2, dist_target = ''; end

close all;
warning('off','MATLAB:unknownObjectIEEE');

ROOT = fileparts(fileparts(mfilename('fullpath')));
try; cd(ROOT); catch; end
addpath(genpath(fullfile(ROOT,'core')));
addpath(fullfile(ROOT,'binaries'));
set(0,'DefaultFigureVisible','off');

cfg   = loadConfig();
N_MAX = 100;
N_VEC = [5 10 20 30 50 75 100];
PARAMS    = {'freq','zS','svp'};
PARAM_LBL = {'Frequency','Source depth zS','SVP shift'};
freq0 = cfg.nominal.freq_Hz;
zS0   = cfg.nominal.zS_m;
geo   = cfg.nominal.geo(:)';

if isempty(gcp('nocreate'))
    try; parpool('local', min(10, maxNumCompThreads)); catch; end
end

for si = 1:numel(cfg.scenarios)
    sc  = cfg.scenarios(si);
    if ~isempty(sc_target) && ~strcmp(sc.name, sc_target), continue; end
    FOM = getFOM(cfg, sc.name);

    fprintf('\n LHS+MoM vs IID+MLE | %s | FOM=%d dB\n', sc.name, FOM);

    raw_cache = fullfile('Cache', sc.name, 'bellhop_raw');
    if ~exist(raw_cache,'dir'), mkdir(raw_cache); end

    sim_nom = {freq0, sc.maxR_m, zS0, 0, 0, 'summer', sc.bathy_type, geo, FOM};
    [TL_dummy, r_grid, z_grid] = bellhopCached(sim_nom, raw_cache);
    r_km = r_grid/1000;  z_m = z_grid;
    [Nz, Nr] = size(TL_dummy);  clear TL_dummy;
    max_depth = sc.maxDepth_m;

    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);
        if ~isempty(dist_target) && ~strcmp(dist.name, dist_target), continue; end
        fprintf('\n  Distribution: %s\n', dist.name);

        fig_dir = fullfile('Methods','LHS_IID_Conv', sc.name, dist.name);
        if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

        B     = computeVarianceBounds(cfg, dist);
        S_lhs = lhsSample(N_MAX, cfg, dist, B, cfg.MC.rng_seed);
        S_iid = iidSample( N_MAX, cfg, dist, B, cfg.MC.rng_seed + 1000);

        for pi = 1:numel(PARAMS)
            param   = PARAMS{pi};
            out_fig = fullfile(fig_dir, sprintf('conv_%s.png', param));
            if isfile(out_fig)
                fprintf('    [SKIP] %s\n', param); continue;
            end

            TL_lhs = buildTLcube(sc, dist, param, S_lhs, 'lhs', N_MAX, ...
                Nz, Nr, max_depth, freq0, zS0, geo, FOM, cfg, raw_cache);
            TL_iid = buildTLcube(sc, dist, param, S_iid, 'iid', N_MAX, ...
                Nz, Nr, max_depth, freq0, zS0, geo, FOM, cfg, raw_cache);

            nv = numel(N_VEC);
            EX_lhs = zeros(Nz,Nr,nv);  Pd_lhs = zeros(Nz,Nr,nv);
            EX_iid = zeros(Nz,Nr,nv);  Pd_iid = zeros(Nz,Nr,nv);

            for ni = 1:nv
                Ni = N_VEC(ni);

                sub = double(TL_lhs(:,:,1:Ni));
                EX_lhs(:,:,ni) = mean(sub,3);
                pix = reshape(permute(sub,[3 1 2]), Ni, Nz*Nr);
                w   = ones(Ni,1)/Ni;
                [~,~,pv] = ln3Moments(pix, FOM, w);
                Pd_lhs(:,:,ni) = reshape(pv, Nz, Nr);

                sub = double(TL_iid(:,:,1:Ni));
                EX_iid(:,:,ni) = mean(sub,3);
                pix = reshape(permute(sub,[3 1 2]), Ni, Nz*Nr);
                pv  = ln3MLE(pix, FOM);
                Pd_iid(:,:,ni) = reshape(pv, Nz, Nr);
            end

            rmse_EX = zeros(1,nv);
            rmse_Pd = zeros(1,nv);
            for ni = 1:nv
                rmse_EX(ni) = sqrt(mean((EX_lhs(:,:,ni)-EX_iid(:,:,ni)).^2,'all'));
                rmse_Pd(ni) = sqrt(mean((Pd_lhs(:,:,ni)-Pd_iid(:,:,ni)).^2,'all'));
            end

            mEX_lhs = squeeze(mean(EX_lhs,[1 2]));
            mEX_iid = squeeze(mean(EX_iid,[1 2]));
            mPd_lhs = squeeze(mean(Pd_lhs,[1 2]));
            mPd_iid = squeeze(mean(Pd_iid,[1 2]));

            fig = figure('Position',[50 50 1100 800]);

            ax1 = subplot(2,3,1);
            plot(ax1,N_VEC,mEX_lhs,'-ob','LineWidth',1.5,'MarkerSize',5,...
                 'DisplayName','LHS + MoM'); hold(ax1,'on');
            plot(ax1,N_VEC,mEX_iid,'-sr','LineWidth',1.5,'MarkerSize',5,...
                 'DisplayName','IID + MLE');
            xlabel(ax1,'N'); ylabel(ax1,'Mean E[TL] (dB)');
            title(ax1,['E[TL] vs N | ' PARAM_LBL{pi}],'Interpreter','none');
            legend(ax1,'Location','best','FontSize',8); grid(ax1,'on');

            ax2 = subplot(2,3,2);
            plot(ax2,N_VEC,mPd_lhs,'-ob','LineWidth',1.5,'MarkerSize',5,...
                 'DisplayName','LHS + MoM'); hold(ax2,'on');
            plot(ax2,N_VEC,mPd_iid,'-sr','LineWidth',1.5,'MarkerSize',5,...
                 'DisplayName','IID + MLE');
            xlabel(ax2,'N'); ylabel(ax2,'Mean P(detect)');
            title(ax2,'P(detect) vs N','Interpreter','none');
            legend(ax2,'Location','best','FontSize',8); grid(ax2,'on');

            ax3 = subplot(2,3,3);
            yyaxis(ax3,'left');
            semilogy(ax3,N_VEC,max(rmse_EX,1e-6),'-ok','LineWidth',1.5,'MarkerSize',5);
            ylabel(ax3,'RMSE E[TL] (dB)');
            yyaxis(ax3,'right');
            semilogy(ax3,N_VEC,max(rmse_Pd,1e-6),'--^k','LineWidth',1.5,'MarkerSize',5);
            ylabel(ax3,'RMSE P(detect)');
            xlabel(ax3,'N');
            title(ax3,'RMSE(LHS vs IID) converges to 0','Interpreter','none');
            grid(ax3,'on');

            ax4 = subplot(2,3,4);
            imagesc(ax4,r_km,z_m,EX_lhs(:,:,end));
            set(ax4,'YDir','reverse'); colorbar(ax4); clim(ax4,[50 150]);
            title(ax4,sprintf('LHS+MoM  E[TL], N=%d',N_MAX),'Interpreter','none');
            xlabel(ax4,'Range (km)'); ylabel(ax4,'Depth (m)');

            ax5 = subplot(2,3,5);
            imagesc(ax5,r_km,z_m,EX_iid(:,:,end));
            set(ax5,'YDir','reverse'); colorbar(ax5); clim(ax5,[50 150]);
            title(ax5,sprintf('IID+MLE  E[TL], N=%d  (FINAL)',N_MAX),'Interpreter','none');
            xlabel(ax5,'Range (km)'); ylabel(ax5,'Depth (m)');

            ax6 = subplot(2,3,6);
            imagesc(ax6,r_km,z_m,EX_lhs(:,:,end)-EX_iid(:,:,end));
            set(ax6,'YDir','reverse'); colorbar(ax6);
            title(ax6,sprintf('Diff LHS-IID (N=%d)',N_MAX),'Interpreter','none');
            xlabel(ax6,'Range (km)'); ylabel(ax6,'Depth (m)');
            colormap(ax6,redblue(256));

            sgtitle(sprintf('LHS+MoM vs IID+MLE | %s | %s | %s', ...
                sc.name, dist.name, PARAM_LBL{pi}), ...
                'Interpreter','none','FontSize',11);

            saveFigPNG(fig, strrep(out_fig,'.png',''));
            drawnow; close(fig);
            fprintf('    Saved: %s\n', out_fig);
        end
    end
end
fprintf('\nrun_lhs_iid_convergence complete\n');
end


function TL_all = buildTLcube(sc, dist, param, S, tag, N, ...
                               Nz, Nr, max_depth, freq0, zS0, geo, FOM, cfg, raw_cache)
    cdir = fullfile('Cache', sc.name, [dist.name '_' tag]);
    if ~exist(cdir,'dir'), mkdir(cdir); end
    tl_file = fullfile(cdir, sprintf('TL_%s.mat', param));

    if isfile(tl_file)
        tmp = load(tl_file,'TL_save','N_saved');
        if tmp.N_saved >= N
            TL_all = double(tmp.TL_save(:,:,1:N));
            fprintf('    [CACHE] %s/%s/%s\n', tag, dist.name, param);
            return;
        end
    end

    do_freq = strcmp(param,'freq');
    do_zS   = strcmp(param,'zS');
    do_svp  = strcmp(param,'svp');

    fprintf('    Running N=%d [%s] %s/%s...\n', N, tag, dist.name, param);
    t0 = tic;  TL_col = cell(N,1);

    parfor i = 1:N
        fi = freq0; if do_freq, fi = S.freq(i); end
        zi = zS0;   if do_zS,   zi = S.zS(i);  end
        sv = 0;     if do_svp,  sv = S.svp(i); end
        sim = {fi, sc.maxR_m, zi, 0, 0, 'summer', sc.bathy_type, geo, FOM};
        tmp_dir = [tempname '_blhp']; mkdir(tmp_dir); prev = cd(tmp_dir);
        if do_svp
            svp_i = makeSVPNoise(sv, max_depth, cfg);
            sim{6} = 'custom';
            TL_col{i} = bellhopCached(sim, raw_cache, 'CustomSVP', svp_i);
        else
            TL_col{i} = bellhopCached(sim, raw_cache);
        end
        cd(prev); rmdir(tmp_dir,'s');
    end

    TL_all = zeros(Nz,Nr,N);
    for i = 1:N, TL_all(:,:,i) = TL_col{i}; end
    fprintf('    Done in %.0fs\n', toc(t0));
    TL_save = single(TL_all);  N_saved = int32(N);
    save(tl_file,'TL_save','N_saved','-v7.3');
end
