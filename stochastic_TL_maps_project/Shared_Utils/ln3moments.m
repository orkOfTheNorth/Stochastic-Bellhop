function [EX, Var, prob] = ln3moments(TL_mat, FOM, w)
% Vectorized 3-parameter lognormal moments for all pixels simultaneously.
%
% Supports optional IS weights for recycled MC distributions.
%
% Inputs:
%   TL_mat — [N × Npix] TL samples
%   FOM    — scalar threshold (dB); needed only when prob output requested
%   w      — [N × 1] normalized IS weights (optional; default = uniform 1/N)
%
% Outputs:
%   EX   — [1 × Npix]  E[TL] under LN3
%   Var  — [1 × Npix]  Var[TL] under LN3
%   prob — [1 × Npix]  P(TL < FOM) = P(detection); NaN when FOM not supplied

TL_mat = double(TL_mat);
N = size(TL_mat, 1);

%% Weights
if nargin < 3 || isempty(w)
    w = ones(N,1) / N;
else
    w = w(:) / sum(w);   % ensure normalized
end
w_col = w;               % [N × 1]

%% Shift parameter: 0.95 × weighted min per pixel
% Weighted quantile approximation: use weighted_min ≈ min of samples with w > 1/(10N)
gam = 0.95 * min(TL_mat, [], 1);   % [1 × Npix]

sh  = max(TL_mat - gam, 1e-6);    % [N × Npix] shifted
lsh = log(sh);

%% Weighted mean and variance of log-shifted values
mu   = w_col' * lsh;                           % [1 × Npix]
res  = lsh - mu;                               % [N × Npix]
sig2 = w_col' * (res .^ 2);                   % [1 × Npix] weighted variance
% Unbias correction for reliability weights: divide by (1 - sum(w²))
sig2 = sig2 / max(1 - sum(w.^2), eps);
sig2 = max(sig2, 1e-12);
sig  = sqrt(sig2);

EX  = gam + exp(mu + sig2/2);
Var = exp(2*mu + sig2) .* (exp(sig2) - 1);

if nargout >= 3
    if nargin < 2 || isempty(FOM)
        prob = nan(size(EX));
    else
        prob             = logncdf(FOM - gam, mu, sig);   % P(TL < FOM) = P(detect)
        prob(FOM <= gam) = 0.0;
    end
end
end
