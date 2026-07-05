function run_ln3_scatter(sc_target, dist_target)
%% run_ln3_scatter  MOM vs MLE P_shadow scatter — fully vectorised.
%
% For every pixel, computes P_shadow via:
%   MLE: gamma=0.95*min; mu=mean(log(Y)); sig=std(log(Y)); Pd=logncdf(FOM-gamma,mu,sig)
%   MOM: same gamma; M1=mean(Y); sig=sqrt(log(M2/M1^2)); mu=log(M1)-sig^2/2
%
% Output: Methods/LN3/<sc>/<dist>/figures/LN3_MOM_vs_MLE_scatter_<param>.png

if nargin < 1, sc_target  = 'deep_water';   end
if nargin < 2, dist_target = 'Normal_5pct'; end

close all;
ROOT = fileparts(fileparts(mfilename('fullpath')));
try; cd(ROOT); catch; end
addpath(genpath(fullfile(ROOT,'core')));
set(0,'DefaultFigureVisible','off');

cfg    = loadConfig();
PARAMS = {'freq','zS','svp'};

for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);
    if ~strcmp(sc.name, sc_target), continue; end
    FOM = getFOM(cfg, sc.name);

    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);
        if ~strcmp(dist.name, dist_target), continue; end

        fig_dir = fullfile('Methods','LN3', sc.name, dist.name, 'figures');
        if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

        for pi = 1:numel(PARAMS)
            param   = PARAMS{pi};
            out_png = fullfile(fig_dir, sprintf('LN3_MOM_vs_MLE_scatter_%s.png', param));
            if isfile(out_png)
                fprintf('[SKIP] %s\n', out_png); continue;
            end

            cache_file = fullfile('Cache', sc.name, dist.name, ...
                                  sprintf('TL_%s.mat', param));
            if ~isfile(cache_file)
                fprintf('[MISS] %s\n', cache_file); continue;
            end

            fprintf('Loading %s...\n', cache_file);
            d    = load(cache_file, 'TL_save');
            TL   = double(d.TL_save);   % [Nz x Nr x N]
            clear d;

            % ── gamma (shared by MLE and MOM, same as run_ln3.m) ─────────
            gam = 0.95 * min(TL, [], 3);          % [Nz x Nr]
            sh  = max(TL - gam, 1e-6);            % [Nz x Nr x N]

            % ── MLE ──────────────────────────────────────────────────────
            lsh     = log(sh);
            mu_mle  = mean(lsh, 3);
            sig_mle = max(std(lsh, 0, 3), 1e-6);
            fom_sh  = max(FOM - gam, 1e-9);
            pd_mle  = normcdf((log(fom_sh) - mu_mle) ./ sig_mle);  % [Nz x Nr]

            % ── MOM ──────────────────────────────────────────────────────
            M1     = mean(sh,   3);
            M2     = mean(sh.^2, 3);
            ratio  = M2 ./ max(M1.^2, 1e-12);
            ratio  = max(ratio, 1+1e-9);           % enforce ratio > 1
            sig_mom = sqrt(log(ratio));
            mu_mom  = log(max(M1,1e-9)) - 0.5 * sig_mom.^2;
            sig_mom = max(sig_mom, 1e-6);
            pd_mom  = normcdf((log(fom_sh) - mu_mom) ./ sig_mom);

            clear TL sh lsh;

            % ── scatter ──────────────────────────────────────────────────
            x = pd_mle(:);
            y = pd_mom(:);

            fig = figure('Position',[50 50 500 480]);
            scatter(x, y, 3, 'filled', 'MarkerFaceAlpha',0.15, ...
                    'MarkerFaceColor',[0.2 0.4 0.8]);
            hold on;
            plot([0 1],[0 1],'k--','LineWidth',1.2);
            xlabel('P_{shadow} — MLE');
            ylabel('P_{shadow} — MOM');
            title(sprintf('LN3 MOM vs MLE | %s | %s | %s', ...
                          sc.name, dist.name, param), 'Interpreter','none');
            axis equal; axis([0 1 0 1]); grid on;
            rmse_val = sqrt(mean((x - y).^2));
            text(0.05, 0.93, sprintf('RMSE = %.4f', rmse_val), ...
                 'FontSize',9,'Color',[0.5 0 0]);

            print(fig, out_png, '-dpng', '-r150');
            close(fig);
            fprintf('Saved: %s\n', out_png);
        end
    end
end
fprintf('run_ln3_scatter done.\n');
end
