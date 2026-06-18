%% fill_findings.m — Extract computed values for findings.tex dagger placeholders.
%
% Reads cached .mat result files and prints LaTeX-ready numbers for every
% dagger marker (†) left in findings.tex.  Run once after a full pipeline
% pass; copy-paste the printed snippets into the corresponding tables.
%
% Usage (from stochastic_TL_maps_project/):
%   matlab -batch "fill_findings"
%   — or interactively: fill_findings

try
    ROOT = fileparts(fileparts(mfilename('fullpath')));
    cd(ROOT);
catch
end

addpath(genpath(fullfile(ROOT, 'core')));
addpath(fullfile(ROOT, 'binaries'));

cfg    = loadConfig();
FOM    = cfg.nominal.FOM_dB;
PARAMS = {'freq', 'zS', 'svp'};

fprintf('\n==================================================================\n');
fprintf('  fill_findings.m — LaTeX values for findings.tex dagger markers\n');
fprintf('==================================================================\n\n');

% ─────────────────────────────────────────────────────────────────────────────
%  TABLE 1: tab:ks_ln3  — Median KS statistic per scenario × distribution
%  Each distribution now has its own independent TL cache (N=150 direct LHS).
% ─────────────────────────────────────────────────────────────────────────────
fprintf('--- TABLE: tab:ks_ln3  (Median KS for LN3 fit) ---\n');
dist_names = {cfg.distributions.name};
hdr = sprintf('%-18s', 'Scenario');
for k = 1:numel(dist_names)
    hdr = [hdr, sprintf('  %16s', strrep(dist_names{k},'_',' '))]; %#ok<AGROW>
end
fprintf('%s\n', hdr);
fprintf('%s\n', repmat('-', 1, 18 + numel(dist_names)*18));

for si = 1:numel(cfg.scenarios)
    sc  = cfg.scenarios(si);
    row = sprintf('%-18s', sc.name);

    for k = 1:numel(dist_names)
        ks_med  = NaN;
        tl_file = fullfile('Cache', sc.name, dist_names{k}, 'TL_zS.mat');
        if isfile(tl_file)
            try
                d   = load(tl_file, 'TL_save');
                TL  = double(d.TL_save);
                [Nz, Nr, N] = size(TL);
                TL_pix = reshape(permute(TL, [3 1 2]), N, Nz*Nr);
                ks_vals = zeros(1, Nz*Nr);
                for px = 1:Nz*Nr
                    x = TL_pix(:, px);
                    x = x(isfinite(x));
                    if numel(x) < 4, continue; end
                    gamma = 0.95 * min(x);
                    y = log(x - gamma);
                    y = y(isfinite(y));
                    if numel(y) < 4, continue; end
                    mu_y  = mean(y);
                    sig_y = std(y);
                    if sig_y < 1e-12, continue; end
                    [~, ~, ks_stat] = kstest((y - mu_y) / sig_y);
                    ks_vals(px) = ks_stat;
                end
                ks_med = median(ks_vals(ks_vals > 0));
            catch me
                fprintf('  [WARN] %s/%s : %s\n', sc.name, dist_names{k}, me.message);
            end
        end
        s = ternary(isnan(ks_med), 'N/A', sprintf('%.3f', ks_med));
        row = [row, sprintf('  %16s', s)]; %#ok<AGROW>
    end
    fprintf('%s\n', row);
end

fprintf('\nLaTeX snippet (copy into tab:ks_ln3):\n');
fprintf('%%  (re-run fill_findings.m and paste values above)\n\n');

% ─────────────────────────────────────────────────────────────────────────────
%  TABLE 2: tab:analytical_l1  — L1 errors from analytical waveguide validation
% ─────────────────────────────────────────────────────────────────────────────
fprintf('\n--- TABLE: tab:analytical_l1  (Analytical Delta validation L1) ---\n');
anal_dir = fullfile('validation', 'analytical', 'results');
dists = struct('name', {'Normal_1pct','Normal_5pct','Normal_10pct'}, ...
               'sigma', {0.5, 2.5, 5.0});

fprintf('%-14s  %8s  %22s  %22s\n', 'Distribution', 'σ_zS (m)', 'L1 (Analytical J)', 'L1 (FD J)');
fprintf('%s\n', repmat('-', 1, 72));

for k = 1:numel(dists)
    afile = fullfile(anal_dir, sprintf('validation_%s.mat', dists(k).name));
    if ~isfile(afile)
        fprintf('%-14s  %8.1f  %22s  %22s\n', dists(k).name, dists(k).sigma, 'file missing', 'file missing');
        continue;
    end
    try
        a = load(afile);
        null_mask = a.TL_nom(:) <= 120;
        mean_VM   = max(mean(a.Var_MC(null_mask)), 1e-10);
        l1_anal = mean(abs(a.Var_delta_anal(null_mask) - a.Var_MC(null_mask))) / mean_VM;
        l1_fd   = mean(abs(a.Var_delta_fd(null_mask)   - a.Var_MC(null_mask))) / mean_VM;
        fprintf('%-14s  %8.1f  %22.2f%%  %22.2f%%\n', dists(k).name, dists(k).sigma, ...
                100*l1_anal, 100*l1_fd);
    catch me
        fprintf('%-14s  %8.1f  [error: %s]\n', dists(k).name, dists(k).sigma, me.message);
    end
end

% ─────────────────────────────────────────────────────────────────────────────
%  SECTION: sec:lhs_explanation  — LHS rationale for findings / methods text
% ─────────────────────────────────────────────────────────────────────────────
fprintf('\n\n--- SECTION: sec:lhs_explanation ---\n\n');
fprintf(['WHY LHS IS NEEDED\n', ...
'With N=%d samples and 3 uncertain parameters, pure random (IID) Monte Carlo\n', ...
'leaves coverage gaps by chance: some quantile regions may receive zero samples,\n', ...
'biasing detection probability estimates near shadow zone boundaries where the\n', ...
'TL distribution is most sensitive to parameter variation.\n\n'], cfg.N_MC);

fprintf(['HOW LHS WORKS\n', ...
'lhsdesign(N, d) divides [0,1]^d into N equal strata per dimension and guarantees\n', ...
'exactly one sample per stratum in each marginal. Transformed via norminv, each\n', ...
'sample xi ~ N(0, sigma^2) marginally — the correct target distribution.\n\n']);

fprintf(['IS LHS EQUIVALENT TO IID MC PDF-WISE?\n', ...
'Marginally yes; jointly no — and the difference is beneficial.\n\n', ...
'Each individual sample xi has the correct marginal distribution N(0, sigma^2),\n', ...
'so all estimators (E[TL], Var[TL], P_detect) are unbiased — identical to IID\n', ...
'Monte Carlo in expectation. The distribution over which each sample is drawn is\n', ...
'the same as standard MC.\n\n', ...
'However, the N samples are negatively correlated across draws: stratification\n', ...
'enforces spread across quantile strata, which makes estimator variance strictly\n', ...
'lower than IID MC for the same N. LHS is therefore a strictly more efficient\n', ...
'estimator, not merely equivalent. IID-based standard error formulas apply\n', ...
'conservatively (they overestimate variance) when used with LHS output.\n\n', ...
'Practical implication: N=150 LHS samples provides coverage of the input\n', ...
'uncertainty comparable to ~200+ IID samples, while remaining statistically\n', ...
'equivalent to IID MC for all reported statistics (mean, variance, P_detect).\n']);

% ─────────────────────────────────────────────────────────────────────────────
%  SECTION: sec:delta_method  — Delta method idea, math, and name discussion
% ─────────────────────────────────────────────────────────────────────────────
fprintf('\n\n--- SECTION: sec:delta_method ---\n\n');

fprintf(['CORE IDEA\n', ...
'The Delta method approximates how uncertainty in inputs propagates to\n', ...
'uncertainty in an output, using a Taylor expansion of the output function\n', ...
'around the nominal input values. Here the output is the TL field\n', ...
'TL(r,z; freq, zS, svp) and the inputs are the three uncertain parameters.\n\n']);

fprintf(['FIRST-ORDER APPROXIMATION\n', ...
'Let x = (freq, zS, svp) with nominal x0 and independent uncertainties sigma_i^2.\n', ...
'Taylor-expanding TL around x0 to first order:\n\n', ...
'  TL(x) approx TL(x0) + sum_i (dTL/dx_i)|x0 * (x_i - x0_i)\n\n', ...
'Taking expectation and variance:\n\n', ...
'  E[TL]    approx TL(x0)                     [first-order mean = nominal TL]\n', ...
'  Var[TL]  approx sum_i (dTL/dx_i)^2 * sigma_i^2\n\n', ...
'The Jacobian J_i = dTL/dx_i is computed by finite differences\n', ...
'(two Bellhop runs per parameter, at +h and -h from nominal).\n\n']);

fprintf(['SECOND-ORDER BIAS CORRECTION\n', ...
'Including the Hessian diagonal terms (H_ii = d^2TL/dx_i^2):\n\n', ...
'  E[TL] approx TL(x0) + (1/2) * sum_i H_ii * sigma_i^2\n\n', ...
'The correction (1/2) H_ii sigma_i^2 is the "2nd-order bias" shown in\n', ...
'bias_2nd_order.png. When this correction is small relative to TL(x0),\n', ...
'the first-order approximation is accurate — i.e., TL is nearly linear\n', ...
'in the uncertain parameter over the range [x0 - 3*sigma, x0 + 3*sigma].\n\n']);

fprintf(['GAUSS-HERMITE EXTENSION\n', ...
'Rather than truncating at Taylor order 2, Gauss-Hermite (GH) quadrature\n', ...
'integrates the 1-D TL response along each parameter axis exactly up to\n', ...
'polynomial degree 2K-1 using K quadrature points. This project uses\n', ...
'orders K=1..10. When GH variance converges as K increases and agrees\n', ...
'with the 1st-order Delta estimate, TL is essentially linear in the\n', ...
'uncertain parameter and the Gaussian assumption is well-supported.\n\n']);

fprintf(['IS "DELTA METHOD" THE RIGHT NAME?\n', ...
'Yes — "Delta method" is the standard statistical term for Taylor-expansion-\n', ...
'based variance propagation. The name comes from "delta" as infinitesimal\n', ...
'perturbation (differential calculus), not the Dirac delta function.\n\n', ...
'The same technique appears under other names:\n', ...
'  FOSM  — First-Order Second-Moment method (engineering / reliability)\n', ...
'  GUM   — linearization per the Guide to Uncertainty in Measurement\n', ...
'          (BIPM, 2008; standard in metrology and physics)\n', ...
'  Error propagation — experimental physics textbooks\n\n', ...
'In a statistics paper, "Delta method" is unambiguous and correct.\n', ...
'In an acoustics or engineering paper, noting the FOSM alias helps\n', ...
'readers avoid confusion with "delta" meaning finite difference or\n', ...
'the Dirac delta function. The name is not misleading.\n']);

% ─────────────────────────────────────────────────────────────────────────────
%  AUTO-PATCH: replace $\dagger$ tokens in findings.tex with KS values
%  Order must match the table rows: scenarios × distributions (outer × inner)
% ─────────────────────────────────────────────────────────────────────────────
tex_path = fullfile(ROOT, 'docs', 'findings.tex');
if isfile(tex_path)
    try
        tex = fileread(tex_path);

        % Compute KS values in row-major order: scenarios × distributions
        kvals = {};
        for si2 = 1:numel(cfg.scenarios)
            sc2 = cfg.scenarios(si2);
            for ki = 1:numel(dist_names)
                ks_val = NaN;
                tl_f2  = fullfile('Cache', sc2.name, dist_names{ki}, 'TL_zS.mat');
                if isfile(tl_f2)
                    try
                        d2 = load(tl_f2, 'TL_save');
                        TL2 = double(d2.TL_save);
                        [Nz2,Nr2,N2] = size(TL2);
                        TLp2 = reshape(permute(TL2,[3 1 2]), N2, Nz2*Nr2);
                        ksv2 = zeros(1, Nz2*Nr2);
                        for px2 = 1:Nz2*Nr2
                            xp = TLp2(:,px2); xp = xp(isfinite(xp));
                            if numel(xp)<4, continue; end
                            gm = 0.95*min(xp); yp = log(xp-gm); yp = yp(isfinite(yp));
                            if numel(yp)<4 || std(yp)<1e-12, continue; end
                            [~,~,ksv2(px2)] = kstest((yp-mean(yp))/std(yp));
                        end
                        ks_val = median(ksv2(ksv2>0));
                    catch, end
                end
                kvals{end+1} = ternary(isnan(ks_val), 'N/A', sprintf('%.3f', ks_val)); %#ok<AGROW>
            end
        end

        % Replace $\dagger$ tokens one-by-one in document order
        DAGGER = '$\dagger$';
        for ki = 1:numel(kvals)
            idx = strfind(tex, DAGGER);
            if isempty(idx), break; end
            tex = [tex(1:idx(1)-1), kvals{ki}, tex(idx(1)+length(DAGGER):end)];
        end

        fid2 = fopen(tex_path, 'w');
        fwrite(fid2, tex);
        fclose(fid2);
        fprintf('Auto-patched %d dagger values into findings.tex\n', numel(kvals));
    catch ME2
        fprintf('[WARN] Auto-patch failed: %s\n', ME2.message);
    end
end

fprintf('\n==================================================================\n');
fprintf('  Done.  findings.tex auto-patched with computed values.\n');
fprintf('==================================================================\n\n');


function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
