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

FOM = 100;   % dB

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
    shd_delta_cheb = D.Cheb_lb >= 0.95;   % Delta Chebyshev 95%
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
saveas(fig_EX, fullfile('figures','compare_EX_diff_all.png'));

figure(fig_Var); subplot(2,4,8); axis off;
sgtitle('MC − Delta: Var[TL] Difference Maps for All 7 Subsets  (red=MC higher, blue=Delta higher)');
saveas(fig_Var, fullfile('figures','compare_Var_diff_all.png'));

figure(fig_Shd); subplot(2,4,8); axis off;
sgtitle(sprintf(['Shadow-Zone Comparison (P(TL>%ddB)≥95%%) | All 7 Subsets\n' ...
    'Blue=Delta-Cheb only | Red=MC-Empirical only | Green=Both agree'],FOM));
saveas(fig_Shd, fullfile('figures','compare_shadow_all.png'));

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

shd_delta = D.Cheb_lb >= 0.95;
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
          mtitles{k}, D.freq0, D.zS0, FOM));
end
sgtitle(sprintf(['Shadow Zone (P(TL>%ddB)≥95%%) – 4 Methods Compared | Full 3-Parameter Perturbation\n' ...
    'Green=shadow (guaranteed), Grey=not shadow'], FOM));
saveas(fig_full, fullfile('figures','compare_shadow_4methods_full.png'));

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
    'f=%dHz±1%%  zS=%.1fm±1%%  Temp±1°C  |  N=%d MC runs'], D.freq0, D.zS0, M.N));

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

saveas(fig_diff, fullfile('figures','compare_EX_Var_full_interactive.png'));

fprintf('\n=== Comparison complete. Results in results/  Figures in figures/ ===\n');

%% ── LOCAL FUNCTIONS ──────────────────────────────────────────────────────────
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
