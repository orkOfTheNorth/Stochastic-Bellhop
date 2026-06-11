function [EX, Var, prob] = ln3moments(TL_mat, FOM)
% Vectorized 3-parameter lognormal moments for all pixels simultaneously.
%
% Replaces the per-pixel ln3fit loop. Uses MLE (N-denominator) to match
% lognfit behaviour. ~400x faster than a scalar loop on 440K pixels.
%
% Inputs:
%   TL_mat — [N × Npix] TL samples (any numeric, cast to double internally)
%   FOM    — scalar threshold (dB); required only when prob output is requested
%
% Outputs:
%   EX   — [1 × Npix]  E[TL] under LN3
%   Var  — [1 × Npix]  Var[TL] under LN3
%   prob — [1 × Npix]  P(TL > FOM); NaN when FOM not supplied

TL_mat = double(TL_mat);
gam  = 0.95 * min(TL_mat, [], 1);            % [1 × Npix] shift
sh   = max(TL_mat - gam, 1e-6);              % [N × Npix] shifted, always > 0
lsh  = log(sh);
mu   = mean(lsh, 1);                          % [1 × Npix] MLE log-mean
sig2 = mean((lsh - mu).^2, 1);               % [1 × Npix] MLE log-variance (N denom)
sig  = sqrt(sig2);

EX  = gam + exp(mu + sig2/2);
Var = exp(2*mu + sig2) .* (exp(sig2) - 1);

if nargout >= 3
    if nargin < 2
        prob = nan(size(EX));
    else
        prob                = 1 - logncdf(FOM - gam, mu, sig);
        prob(FOM <= gam)    = 1.0;
    end
end
end
