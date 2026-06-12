function [prob, gam, mu_ln, sig_ln] = ln3fit(samps, FOM)
% Fit a 3-parameter lognormal to TL samples and return P(TL < FOM) = P(detect).
%
%   gamma   = 0.95 * min(samples)           [location / shift parameter]
%   shifted = samples - gamma               [always > 0]
%   [mu_ln, sig_ln] = lognfit(shifted)      [2-param lognormal on shifted]
%   prob = logncdf(FOM - gamma, mu_ln, sig_ln)   [P(TL < FOM) = P(detect)]
%
% Inputs:
%   samps — (N x 1) vector of TL realisations at one (z, r) pixel
%   FOM   — scalar threshold in dB
%
% Outputs:
%   prob   — P(TL < FOM) = P(detection) in [0, 1]
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
    prob = logncdf(FOM - gam, mu_ln, sig_ln);   % P(TL < FOM) = P(detect)
else
    prob = 0.0;   % FOM below shift → detection impossible
end
end
