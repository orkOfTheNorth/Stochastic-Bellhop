function [prob, mu_ln, sig_ln, gam] = deltaLN3prob(EX, Var, FOM)
% Moment-matched LN3 P(shadow) from Delta method moments — no MC needed.
%
% Given E[TL] and Var[TL] (from Delta), fits a 3-parameter lognormal by
% matching the first two moments, using a physics-motivated shift parameter:
%
%   γ = max(0, E[TL] − 4·σ[TL])       (4-sigma lower bound on TL)
%
% Moment-matching on the shifted variable X' = TL − γ  ~  LN(μ,σ²):
%
%   CV²   = Var / (E[TL] − γ)²
%   σ_ln  = sqrt(log(1 + CV²))
%   μ_ln  = log(E[TL] − γ) − σ_ln²/2
%   P(shadow) = 1 − logncdf(FOM − γ, μ_ln, σ_ln)
%
% Only pixels where E[TL] > FOM get a non-zero probability — consistent
% with chebyshevBound convention (detection zone → prob = 0).

gam    = max(0, EX - 4.*sqrt(max(Var, 0)));
mu_eff = max(EX - gam, 1e-6);

sig_ln = sqrt(log(1 + Var ./ mu_eff.^2));
mu_ln  = log(mu_eff) - sig_ln.^2 / 2;

prob      = zeros(size(EX));
shadow    = (EX - FOM) > 0;                  % detection zone → stays 0
fom_shift = FOM - gam(shadow);
mu_s      = mu_ln(shadow);
sig_s     = sig_ln(shadow);

p       = ones(sum(shadow(:)), 1);           % default: FOM <= gamma → P=1
above   = fom_shift > 0;
p(above) = 1 - logncdf(fom_shift(above), mu_s(above), sig_s(above));

prob(shadow) = min(max(p, 0), 1);
end
