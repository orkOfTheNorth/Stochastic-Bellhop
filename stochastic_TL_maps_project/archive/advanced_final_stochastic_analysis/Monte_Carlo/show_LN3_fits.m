%% show_LN3_fits.m  –  LN3 fit quality diagnostic  (2026-05-26)
%
% Loads the full 3-parameter MC TL sample cube and generates:
%   1. LN3_grid_fits  – 4×5 tiled histogram + LN3 PDF at 20 sample points
%   2. LN3_KS_map     – spatial heatmap of KS statistic (fit quality map)
%
% KS statistic = max|F_empirical − F_LN3|  (lower = better fit)
% Saved to:  Monte_Carlo/LN3_fits/

clear; close all; clc;
try, cd(fileparts(mfilename('fullpath'))); catch; end

out_dir = 'LN3_fits';
if ~exist(out_dir,'dir'), mkdir(out_dir); end

FOM = 100;

%% ── LOAD TL SAMPLES ─────────────────────────────────────────────────────────
tl_file = fullfile('TL_cache','TL_zS_freq_temp.mat');
if isfile(tl_file)
    tmp    = load(tl_file);
    TL_all = double(tmp.TL_save);
    N      = tmp.N_saved;
    fprintf('Loaded TL_cache: N=%d samples\n', N);
else
    tmp    = load(fullfile('results','MC_zS_freq_temp.mat'));
    TL_all = tmp.TL_all;
    N      = tmp.N;
    fprintf('Loaded MC results: N=%d samples\n', N);
end

mc   = load(fullfile('results','MC_zS_freq_temp.mat'), 'r_km','z_m','MC_EX');
r_km = mc.r_km;
z_m  = mc.z_m;
[Nz, Nr, ~] = size(TL_all);

%% ── 1. TILED GRID OF FITS (5 rows depth × 5 cols range = 25 points) ─────────
nr_pts = 5;   nz_pts = 5;
ri_pts = round(linspace(ceil(Nr*0.06), Nr-1, nr_pts));
zi_pts = round(linspace(ceil(Nz*0.06), Nz-1, nz_pts));

fig1 = figure('Name','LN3 Fit Quality – Grid Sample', ...
              'Position',[30 30 1600 1200]);
sp = 0;
for zi = zi_pts
    for ri = ri_pts
        sp = sp + 1;
        samps = double(squeeze(TL_all(zi,ri,:)));

        subplot(nz_pts, nr_pts, sp);
        try
            [ln3_p, gam, mu_ln, sig_ln] = ln3prob(samps, FOM);
            emp_p = mean(samps > FOM);
            ks    = ksstat(samps, gam, mu_ln, sig_ln);

            % 95% empirical run interval: range that 95% of MC samples fell in
            p025 = prctile(samps, 2.5);
            p975 = prctile(samps, 97.5);

            % PDF evaluation grid and LN3 mean curve
            x     = linspace(min(samps)*0.97, max(samps)*1.03, 200);
            y_ln3 = lognpdf(x - gam, mu_ln, sig_ln);

            % LN3 model 95% CI via MLE asymptotic standard errors:
            %   se(mu_ln)  = sig_ln / sqrt(N)
            %   se(sig_ln) = sig_ln / sqrt(2*N)
            % Envelope of PDFs at all 4 corners of the (mu, sig) CI box.
            se_mu = sig_ln / sqrt(N);
            se_sg = sig_ln / sqrt(2*N);
            mu_lo = mu_ln - 1.96*se_mu;   mu_hi = mu_ln + 1.96*se_mu;
            sg_lo = max(sig_ln - 1.96*se_sg, 1e-4);
            sg_hi = sig_ln + 1.96*se_sg;
            all_pdfs = [lognpdf(x-gam, mu_lo, sg_lo);
                        lognpdf(x-gam, mu_lo, sg_hi);
                        lognpdf(x-gam, mu_hi, sg_lo);
                        lognpdf(x-gam, mu_hi, sg_hi)];
            ci_lo = min(all_pdfs, [], 1);
            ci_hi = max(all_pdfs, [], 1);

            hold on;
            % ① LN3 model CI band — light blue fill, drawn first (sits behind bars)
            fill([x, fliplr(x)], [ci_lo, fliplr(ci_hi)], ...
                [0.20 0.45 0.85], 'FaceAlpha', 0.22, 'EdgeColor','none');
            % ② Histogram (empirical PDF)
            histogram(samps, 15, 'Normalization','pdf', ...
                'FaceColor',[0.35 0.60 0.88],'EdgeColor','none','FaceAlpha',0.75);
            % ③ LN3 mean fit line
            plot(x, y_ln3, 'b-', 'LineWidth', 1.0);
            % ④ 95% run range: dotted grey verticals at 2.5th and 97.5th pct
            xline(p025, ':', 'Color',[0.40 0.40 0.40], 'LineWidth', 1.0, ...
                'HandleVisibility','off');
            xline(p975, ':', 'Color',[0.40 0.40 0.40], 'LineWidth', 1.0, ...
                'DisplayName','95% of runs');
            % ⑤ FOM threshold
            xline(FOM, 'k--', 'LineWidth', 1.3);
            hold off;
            ylim([0, max(y_ln3) * 1.45]);

            ttl = sprintf('R=%.0fkm  Z=%.0fm\nKS=%.3f  P_{LN3}=%.0f%%  P_{emp}=%.0f%%', ...
                          r_km(ri), z_m(zi), ks, ln3_p*100, emp_p*100);
            if ks < 0.10,     tc = [0.0 0.55 0.0];
            elseif ks < 0.15, tc = [0.8 0.5  0.0];
            else,             tc = [0.8 0.1  0.1];
            end
            title(ttl, 'FontSize', 6.5, 'Color', tc);
        catch
            title(sprintf('R=%.0fkm Z=%.0fm\n(fit failed)', r_km(ri), z_m(zi)), ...
                  'FontSize',6.5,'Color',[0.6 0.6 0.6]);
        end

        set(gca,'FontSize',6);
        if sp == 1
            legend({'LN3 model CI (95%)', 'Empirical PDF', ...
                    'LN3 fit', '95% of runs', sprintf('FOM=%d dB',FOM)}, ...
                   'Location','northwest','FontSize',5.5);
        end
        if sp > (nz_pts-1)*nr_pts, xlabel('TL (dB)','FontSize',7); end
        if mod(sp-1,nr_pts)==0
            ylabel({'Probability density','(prob per dB)'},'FontSize',7);
        end
    end
end
sgtitle(sprintf(['LN3 Fit Quality – Full 3-Param MC | N=%d\n' ...
    'P_{LN3} = shadow probability from fitted LN3 distribution    ' ...
    'P_{emp} = fraction of %d samples above FOM (direct count, no model assumed)\n' ...
    'Y-axis = probability DENSITY (prob per dB), not probability — integrates to 1 over all TL\n' ...
    'Title colour: green=good (KS<0.10)  orange=ok (KS<0.15)  red=poor (KS>=0.15)'], N, N), ...
    'FontSize', 9);
saveFig(fig1, fullfile(out_dir,'LN3_grid_fits'));
close(fig1);
fprintf('Saved LN3_grid_fits\n');

%% ── 2. SPATIAL KS-STAT MAP (every 5th grid point) ───────────────────────────
stride   = 5;
zi_map   = 1:stride:Nz;
ri_map   = 1:stride:Nr;
nzi      = numel(zi_map);
nri      = numel(ri_map);
ks_map_s = zeros(nzi, nri);

fprintf('Computing KS map (%d×%d points) ...\n', nzi, nri);
for ii = 1:nzi
    for jj = 1:nri
        s = double(squeeze(TL_all(zi_map(ii), ri_map(jj), :)));
        try
            [~, gam, mu_ln, sig_ln] = ln3prob(s, FOM);
            ks_map_s(ii,jj)         = ksstat(s, gam, mu_ln, sig_ln);
        catch
            ks_map_s(ii,jj) = NaN;
        end
    end
    if mod(ii,20)==0 || ii==nzi
        fprintf('  row %d/%d\n', ii, nzi);
    end
end

fig2 = figure('Name','LN3 KS-Stat Spatial Map','Position',[60 60 1050 620]);
pcolor(r_km(ri_map), z_m(zi_map), ks_map_s); shading interp;
set(gca,'YDir','reverse');
cmap = flipud(hot(256));
colormap(gca, cmap); cb = colorbar;
clim([0 0.20]);
ylabel(cb, 'KS statistic (lower = better LN3 fit)');
xlabel('Range (km)'); ylabel('Depth (m)');
title(sprintf(['LN3 Fit Quality Map | KS = max|F_{emp}−F_{LN3}| | N=%d  stride=%d pts\n' ...
    'Dark (low KS) = good fit     Light/white (KS≥0.20) = poor fit'], N, stride));
saveFig(fig2, fullfile(out_dir,'LN3_KS_map'));
close(fig2);
fprintf('Saved LN3_KS_map\n');

%% ── SUMMARY ────────────────────────────────────────────────────────────────
ks_v = ks_map_s(~isnan(ks_map_s));
fprintf('\n=== LN3 Fit Summary ===\n');
fprintf('  KS < 0.10 (good):       %5.1f%% of sampled points\n', 100*mean(ks_v < 0.10));
fprintf('  KS < 0.15 (acceptable): %5.1f%% of sampled points\n', 100*mean(ks_v < 0.15));
fprintf('  Mean KS : %.4f\n', mean(ks_v));
fprintf('  Max  KS : %.4f\n', max(ks_v));
fprintf('\nAll figures saved to %s/\n', out_dir);

%% ── LOCAL FUNCTIONS ──────────────────────────────────────────────────────────
function [prob, gam, mu_ln, sig_ln] = ln3prob(samps, FOM)
    samps  = double(samps(:));
    gam    = min(samps) * 0.95;
    sh     = samps - gam;
    sh(sh<=0) = 1e-6;
    params = lognfit(sh);
    mu_ln  = params(1);
    sig_ln = params(2);
    if FOM > gam
        prob = 1 - logncdf(FOM - gam, mu_ln, sig_ln);
    else
        prob = 1.0;
    end
end

function ks = ksstat(samps, gam, mu_ln, sig_ln)
    samps  = sort(double(samps(:)));
    n      = numel(samps);
    F_emp  = (1:n)' / n;
    F_ln3  = logncdf(samps - gam, mu_ln, sig_ln);
    ks     = max(abs(F_emp - F_ln3));
end

function saveFig(fig, base_path)
    print(fig, base_path, '-dpng', '-r300');
    try
        exportgraphics(fig, [base_path '.pdf'], 'ContentType','image', 'Resolution',300);
    catch
    end
end
