function nll = ln3NLL(samps, gam, mu, sig)
% Mean negative log-likelihood per sample of a fitted 3-parameter lognormal.
%
%   X - gam ~ LN2(mu, sig)   =>   log p(x) = -log(x-gam) - log(sig) - 0.5*log(2*pi)
%                                            - (log(x-gam)-mu)^2 / (2*sig^2)
%
% Samples with x <= gam fall outside the model's support (probability 0);
% they are penalized with a large finite value (1e6) rather than -Inf so
% that averages/medians over pixels stay finite and comparable.
%
% Inputs:
%   samps — (N x 1) vector, or (P x N) matrix (rows = pixels)
%   gam, mu, sig — fitted LN3 parameters: scalar (if samps is a vector),
%                  or (P x 1) / (1 x P) to match rows (if samps is a matrix)
%
% Output:
%   nll — mean negative log-likelihood per sample: scalar, or (P x 1) vector

if isvector(samps)
    samps = samps(:)';   % 1 x N
    gam = gam(1); mu = mu(1); sig = sig(1);
end
gam = gam(:);  mu = mu(:);  sig = max(sig(:), 1e-8);   % P x 1

PENALTY = 1e6;
sh = samps - gam;                       % P x N (broadcast)
valid = sh > 0;

ly = nan(size(sh));
ly(valid) = log(sh(valid));

logp = -ly - log(sig) - 0.5*log(2*pi) - (ly - mu).^2 ./ (2*sig.^2);
logp(~valid) = -inf;

nll_terms = -logp;
nll_terms(~isfinite(nll_terms)) = PENALTY;

nll = mean(nll_terms, 2);
if isscalar(nll), nll = nll(1); end
end
