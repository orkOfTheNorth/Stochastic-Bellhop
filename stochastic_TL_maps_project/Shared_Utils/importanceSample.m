function [stats, ESS] = importanceSample(TL_save, x_samples, sigma_large, sigma_target, FOM)
%IMPORTANCESAMPLE  IS-weighted statistics from a broad-distribution MC run.
%
%  TL_save     [Nz x Nr x N]  — TL field samples (from P = N(0,sigma_large))
%  x_samples   [N x 1]        — drawn perturbation values under P
%  sigma_large scalar          — std dev of base distribution P
%  sigma_target scalar         — std dev of target distribution Q
%  FOM         scalar          — detection threshold (dB); P_detect = P(TL < FOM)
%
%  stats struct fields:
%    EX        [Nz x Nr]   — IS-weighted mean TL
%    Var       [Nz x Nr]   — IS-weighted variance of TL
%    P_detect  [Nz x Nr]   — IS-weighted P(TL < FOM)
%    w_norm    [N x 1]     — normalized weights (for external use)
%  ESS   scalar  — Kish effective sample size

    N = numel(x_samples);

    %% Importance weights w_i = Q(x_i)/P(x_i)  (log-space for stability)
    % Q,P both N(0,σ²): log w_i = x_i²/2 · (1/σ_large² - 1/σ_target²)
    log_w = (x_samples(:).^2) / 2 .* (1/sigma_large^2 - 1/sigma_target^2);

    % Clip log weights to avoid Inf (far-tail samples)
    log_w = log_w - max(log_w);         % shift for numerical stability
    w     = exp(log_w);
    w_sum = sum(w);

    if w_sum < eps
        warning('importanceSample:zeroWeights', ...
                'All IS weights are zero — target σ much smaller than base σ. Returning NaN.');
        stats.EX       = nan(size(TL_save,1), size(TL_save,2));
        stats.Var      = nan(size(TL_save,1), size(TL_save,2));
        stats.P_detect = nan(size(TL_save,1), size(TL_save,2));
        stats.w_norm   = zeros(N,1);
        ESS = 0;
        return;
    end

    w_norm = w / w_sum;                 % normalized: sum = 1

    %% Kish effective sample size
    ESS = 1 / sum(w_norm.^2);

    if ESS < 30
        warning('importanceSample:lowESS', ...
                'ESS = %.1f (< 30). IS estimates unreliable for σ_target=%.3f, σ_large=%.3f.', ...
                ESS, sigma_target, sigma_large);
    end

    %% Weighted statistics over the N dimension (dim 3)
    % Reshape weights for broadcasting: [1 x 1 x N]
    w3 = reshape(w_norm, 1, 1, N);

    % Weighted mean
    EX = sum(TL_save .* w3, 3);

    % Weighted variance  (reliability-weighted, unbiaed correction 1/(1-sum(w²)))
    dTL    = TL_save - EX;
    Var    = sum(w3 .* dTL.^2, 3) / (1 - sum(w_norm.^2));

    % Weighted P(TL < FOM) = detection probability
    P_detect = sum(w3 .* (TL_save < FOM), 3);

    stats.EX       = EX;
    stats.Var      = Var;
    stats.P_detect = P_detect;
    stats.w_norm   = w_norm;
end
