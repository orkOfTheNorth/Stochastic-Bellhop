function prob = ln3prob(samps, FOM)
% 3-parameter lognormal P(TL > FOM).
% Threshold γ = 0.95 * min(samples).  Fits lognormal to (TL − γ).
    samps = double(samps(:));
    gam   = min(samps) * 0.95;
    sh    = samps - gam;
    sh(sh <= 0) = 1e-6;
    p = lognfit(sh);
    if FOM > gam
        prob = 1 - logncdf(FOM - gam, p(1), p(2));
    else
        prob = 1.0;
    end
end
