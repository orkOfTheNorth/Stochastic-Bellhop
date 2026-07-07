function ks = ksStatLN3(samps, gam, mu_ln, sig_ln)
% Kolmogorov-Smirnov statistic between a fitted LN3 CDF and the empirical
% CDF of the samples. Extracted from pipeline/run_ln3.m (was a file-local
% function there, uncallable from other scripts such as
% scripts/compare_ln3_mom_vs_mle.m) so it can be shared.
samps  = sort(double(samps(:)));
n      = numel(samps);
F_emp  = (1:n)' / n;
sh     = samps - gam;
sh(sh <= 0) = 1e-9;
F_ln3  = logncdf(sh, mu_ln, sig_ln);
ks     = max(abs(F_emp - F_ln3));
end
