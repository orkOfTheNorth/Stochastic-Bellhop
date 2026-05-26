%% run_MC.m  –  Monte Carlo UQ  (2026-05-26)
%
% Runs N=200 uniform Monte Carlo samples for each of the 7 parameter subsets.
% For each subset, only the active parameters are varied; inactive ones
% remain at their nominal values.
%
% For the FULL 3-parameter subset, additionally stores the complete TL
% sample cube (TL_all) so the interactive viewer (show_MC.m) can compute
% on-demand LN3 fits and histograms at any grid point.
%
% Saved per subset:  results/MC_<subset>.mat
%   MC_EX     – empirical mean TL        [Nz × Nr]
%   MC_Var    – empirical variance TL    [Nz × Nr]
%   MC_PrFOM  – empirical P(TL>FOM)      [Nz × Nr]
%   Cheb_lb   – Cantelli lb P(TL>FOM)   [Nz × Nr]
%   r_grid / z_grid / r_km / z_m / FOM / N
%
% FULL subset only also saves:
%   TL_all    – full sample cube         [Nz × Nr × N]

clear; close all; clc; warning('off');
try, cd(fileparts(mfilename('fullpath'))); catch; end
addpath(genpath('../../uniform source depth analasis 13.5/code/Functions'));
if ~exist('results','dir'), mkdir('results'); end
if ~exist('figures','dir'), mkdir('figures'); end

%% ── PARAMETERS ──────────────────────────────────────────────────────────────
freq0 = 10000;   maxR = 50000;   zS0 = 5;   FOM = 100;
geo   = [0.989*1500, 1.63, 0.07];
N     = 200;     % MC sample count

% Uniform bounds
freq_bnd = [freq0*0.99, freq0*1.01];
zS_bnd   = [zS0*0.99,  zS0*1.01];
temp_bnd = [-1, 1];

%% ── SUBSET DEFINITIONS ───────────────────────────────────────────────────────
% active flags = [freq_active, zS_active, temp_active]
subset_names  = {'zS','freq','temp','zS_freq','zS_temp','freq_temp','zS_freq_temp'};
subset_labels = {'z_S only','Freq only','Temp only', ...
                 'z_S + Freq','z_S + Temp','Freq + Temp','z_S + Freq + Temp'};
active = {[0 1 0],[1 0 0],[0 0 1],[1 1 0],[0 1 1],[1 0 1],[1 1 1]};

%% ── PRE-GENERATE SAMPLES (Latin Hypercube over all 3 dims) ──────────────────
rng(42);  % reproducibility
lhs = lhsdesign(N, 3);   % N×3  in [0,1]
% Map to physical values
samp_freq = freq_bnd(1) + lhs(:,1) * diff(freq_bnd);   % Hz
samp_zS   = zS_bnd(1)   + lhs(:,2) * diff(zS_bnd);    % m
samp_temp = temp_bnd(1)  + lhs(:,3) * diff(temp_bnd); % °C shift

%% ── DUMMY RUN TO GET GRID DIMENSIONS ────────────────────────────────────────
set(0,'DefaultFigureVisible','off');
sim_nom = {freq0, maxR, zS0, 0, 0, "summer", "const_35", geo, FOM};
[TL_dummy, r_grid, z_grid] = simpleBellhopHazat(sim_nom);
close all;
r_km = r_grid / 1000;
z_m  = z_grid;
[Nz, Nr] = size(TL_dummy);

%% ── MAIN LOOP OVER SUBSETS ───────────────────────────────────────────────────
for s = 1:7
    fl  = active{s};
    sn  = subset_names{s};
    fprintf('\n=== Subset %d/7: [%s] ===\n', s, subset_labels{s});

    TL_all = zeros(Nz, Nr, N);

    for i = 1:N
        if mod(i,20)==0 || i==1
            fprintf('  Run %d/%d ...\n', i, N);
        end

        % Frequency (active or nominal)
        fi = freq0;
        if fl(1), fi = samp_freq(i); end

        % Source depth (active or nominal)
        zi = zS0;
        if fl(2), zi = samp_zS(i); end

        % Temperature shift (active or nominal=0)
        ti = 0;
        if fl(3), ti = samp_temp(i); end

        % Build sim_pars
        sim = {fi, maxR, zi, 0, 0, "summer", "const_35", geo, FOM};
        if fl(3)  % custom SVP
            sim{6} = "custom";
            [TL,~,~] = simpleBellhopHazat(sim, 'CustomSVP', makeSVP(ti));
        else
            [TL,~,~] = simpleBellhopHazat(sim);
        end
        TL_all(:,:,i) = TL;
        close all;
    end

    % Statistics
    MC_EX    = mean(TL_all, 3);
    MC_Var   = var(TL_all, 0, 3);
    MC_PrFOM = mean(TL_all > FOM, 3);   % empirical P(TL > FOM)

    % Cantelli lower bound (same formula as Delta method)
    d  = MC_EX - FOM;
    lb = zeros(Nz, Nr);
    mask = d > 0;
    lb(mask) = max(0,  1 - MC_Var(mask) ./ (MC_Var(mask) + d(mask).^2));
    Cheb_lb = lb;

    % Save (full subset: also save TL_all for interactive viewer)
    fname = fullfile('results', sprintf('MC_%s.mat', sn));
    if s == 7
        save(fname, 'MC_EX','MC_Var','MC_PrFOM','Cheb_lb', ...
                    'TL_all','r_grid','z_grid','r_km','z_m','FOM','N', ...
                    'freq_bnd','zS_bnd','temp_bnd','freq0','zS0');
    else
        save(fname, 'MC_EX','MC_Var','MC_PrFOM','Cheb_lb', ...
                    'r_grid','z_grid','r_km','z_m','FOM','N');
    end
    fprintf('  Saved %s\n', fname);

    % Quick summary figure (saved, not interactive)
    fig = figure('Name',sprintf('MC %s',subset_labels{s}),'Position',[50 50 1400 400]);

    subplot(1,3,1);
    pcolor(r_km, z_m, MC_EX); shading interp; set(gca,'YDir','reverse');
    colormap(gca,jet); colorbar; clim([50 150]);
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('MC E[TL] | %s | N=%d, f=%dHz, zS=%.1fm',subset_labels{s},N,freq0,zS0));

    subplot(1,3,2);
    pcolor(r_km, z_m, MC_Var); shading interp; set(gca,'YDir','reverse');
    colormap(gca,hot); colorbar;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('MC Var[TL] | %s',subset_labels{s}));

    subplot(1,3,3);
    pcolor(r_km, z_m, MC_PrFOM*100); shading interp; set(gca,'YDir','reverse');
    colormap(gca,parula); colorbar; clim([0 100]);
    hold on;
    contour(r_km, z_m, MC_PrFOM, [0.95 0.95], 'w--','LineWidth',2);
    hold off;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('Empirical P(TL>%ddB) %% | %s | white=95%%',FOM,subset_labels{s}));

    saveas(fig, fullfile('figures', sprintf('MC_%s_summary.png', sn)));
    close(fig);
end

set(0,'DefaultFigureVisible','on');
fprintf('\n=== Monte Carlo complete. Open show_MC.m for interactive viewer. ===\n');

%% ── LOCAL FUNCTION ───────────────────────────────────────────────────────────
function svp = makeSVP(temp_shift)
    depths  = linspace(0, 5000, 1000);
    base_td = [0 30 180 400 5000; 25 25 17 13.6 13.6]';
    t  = interp1(base_td(:,1), base_td(:,2) + temp_shift, depths);
    S  = 37;
    sv = 1499.2 + 4.6*t - 0.055*t.^2 + 0.00029*t.^3 + ...
         (1.34 - 0.01*t).*(S-35) + 0.016*depths;
    svp = [depths.' sv.'];
end
