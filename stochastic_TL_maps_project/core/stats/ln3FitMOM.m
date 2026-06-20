function [prob, gam, mu_ln, sig_ln] = ln3FitMOM(samps, FOM)
% Fit 3-parameter lognormal via TRUE Method of Moments (MOM).
% Compare with ln3Fit which uses gamma=0.95*min + MLE via lognfit.
%
% MOM equations for LN3 (X-gamma ~ LN(mu, sigma)):
%   gamma = 0.95 * min(samps)   [same heuristic shift — isolates MOM vs MLE]
%   shifted Y = samps - gamma
%   MOM for LN2(mu, sigma):
%     sigma_MOM^2 = log(E[Y^2] / E[Y]^2)
%     mu_MOM = log(E[Y]) - sigma_MOM^2 / 2
%   where E[Y] and E[Y^2] are estimated from sample moments.
%
% Outputs: same signature as ln3Fit for drop-in comparison.

samps = double(samps(:));
gam   = 0.95 * min(samps);
sh    = samps - gam;
sh(sh <= 0) = 1e-6;

M1 = mean(sh);
M2 = mean(sh.^2);

if M2 <= M1^2 || M1 <= 0
    % Degenerate case: fall back to MLE
    p = lognfit(sh);
    mu_ln  = p(1);
    sig_ln = p(2);
else
    sig_ln = sqrt(log(M2 / M1^2));
    mu_ln  = log(M1) - 0.5 * sig_ln^2;
end

if FOM > gam
    prob = logncdf(FOM - gam, mu_ln, sig_ln);
else
    prob = 0.0;
end
end
