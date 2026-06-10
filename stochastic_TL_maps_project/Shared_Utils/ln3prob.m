function prob = ln3prob(samps, FOM)
% 3-parameter lognormal P(TL > FOM) — minimal wrapper around ln3fit.
prob = ln3fit(samps, FOM);
end
