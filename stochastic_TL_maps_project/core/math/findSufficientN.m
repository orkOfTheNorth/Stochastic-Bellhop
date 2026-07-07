function N_suff = findSufficientN(N_VEC, metric_vec, tol)
% Smallest N in N_VEC beyond which the (monotonically decreasing) metric
% stays within a relative tolerance `tol` of its value at N_VEC(end)
% (the largest N available, treated as the "converged" reference).
%
% Answers the practical question: "how many samples/runs do we actually
% need?" — the smallest N past which more samples buy negligible accuracy.
%
% Inputs:
%   N_VEC      — increasing vector of sample sizes, e.g. [10 20 ... 1000]
%   metric_vec — matching error metric per N (same length, decreasing-ish;
%                NaN entries allowed, e.g. below a regression threshold)
%   tol        — relative tolerance, default 0.10 (10% above the N_VEC(end)
%                value counts as "practically converged")
%
% Output:
%   N_suff — smallest N_VEC(i) such that metric_vec(i) <= (1+tol)*ref for
%            all i' >= i up to the end (monotone-from-here check, so a
%            single noisy dip doesn't falsely qualify); NaN if never reached
%            or if metric_vec(end) is NaN.

if nargin < 3 || isempty(tol), tol = 0.10; end

N_VEC = N_VEC(:)';  metric_vec = metric_vec(:)';
valid = isfinite(metric_vec);
if ~any(valid) || ~valid(end)
    N_suff = NaN;
    return;
end

ref = metric_vec(end);
thresh = (1 + tol) * ref;

n = numel(metric_vec);
N_suff = NaN;
for i = 1:n
    if ~valid(i), continue; end
    if all(metric_vec(i:end) <= thresh | ~valid(i:end))
        N_suff = N_VEC(i);
        break;
    end
end
end
