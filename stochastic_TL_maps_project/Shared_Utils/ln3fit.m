function [prob, gam, mu_ln, sig_ln] = ln3fit(samps, FOM)
% Fit a 3-parameter lognormal to TL samples and return P(TL > FOM).
%
%   gamma   = 0.95 * min(samples)           [location / shift parameter]
%   shifted = samples - gamma               [always > 0]
%   [mu_ln, sig_ln] = lognfit(shifted)      [2-param lognormal on shifted]
%   prob = 1 - logncdf(FOM - gamma, mu_ln, sig_ln)
%
% Inputs:
%   samps — (N x 1) vector of TL realisations at one (z, r) pixel
%   FOM   — scalar threshold in dB
%
% Outputs:
%   prob   — P(TL > FOM) in [0, 1]
%   gam    — shift parameter (dB)
%   mu_ln  — lognormal location after shift
%   sig_ln — lognormal scale after shift

samps = double(samps(:));
gam   = 0.95 * min(samps);
sh    = samps - gam;
sh(sh <= 0) = 1e-6;

p = lognfit(sh);        % [mu_ln, sig_ln]
mu_ln  = p(1);
sig_ln = p(2);

if FOM > gam
    prob = 1 - logncdf(FOM - gam, mu_ln, sig_ln);
else
    prob = 1.0;
end
end
