function prob = ln3MLE(X, FOM)
% Fit 3-parameter lognormal to TL samples via MLE; return P(TL < FOM).
%
% X   — [N x Npix] TL samples (dB)
% FOM — scalar detection threshold (dB)
% prob — [1 x Npix] P(TL < FOM)
%
% Gamma found by grid-search profile MLE, fully vectorised across pixels.

[N, Npix] = size(X);
min_X = min(X, [], 1);                       % [1 x Npix]

fracs   = linspace(0.01, 0.995, 80);
best_LL = -inf(1, Npix);
best_gam = 0.95 * min_X;

for f = fracs
    gam = f * min_X;
    Y   = log(X - gam);                      % [N x Npix]
    mu  = mean(Y, 1);
    sig = sqrt(mean((Y - mu).^2, 1));
    sig = max(sig, 1e-8);
    LL  = -N*log(sig) - sum(log(X - gam), 1);
    better = LL > best_LL;
    best_LL(better)  = LL(better);
    best_gam(better) = gam(better);
end

Y_best  = log(X - best_gam);
mu_mle  = mean(Y_best, 1);
sig_mle = max(sqrt(mean((Y_best - mu_mle).^2, 1)), 1e-8);

log_arg = FOM - best_gam;
log_arg(log_arg <= 0) = NaN;
z    = (log(log_arg) - mu_mle) ./ sig_mle;
prob = normcdf(z);
prob(isnan(z)) = 0;
end
