function S = lhsSample(N, cfg, dist, B, seed)
% Generate N LHS samples for [freq, zS, svp_level] given a distribution.
%
% Inputs:
%   N    — number of samples
%   cfg  — struct from loadConfig()
%   dist — distribution config entry (has .type)
%   B    — bounds struct from computeVarianceBounds()
%   seed — RNG seed for reproducibility (default: cfg.MC.rng_seed)
%
% Output S (struct):
%   S.freq   [N x 1] Hz
%   S.zS     [N x 1] m
%   S.svp    [N x 1] degC

if nargin < 5
    seed = cfg.MC.rng_seed;
end
rng(seed);
lhs = lhsdesign(N, 3);   % N x 3, stratified uniform [0,1]

freq0 = cfg.nominal.freq_Hz;
zS0   = cfg.nominal.zS_m;

switch lower(dist.type)
    case 'uniform'
        S.freq = B.freq_bnd(1) + lhs(:,1) * diff(B.freq_bnd);
        S.zS   = B.zS_bnd(1)   + lhs(:,2) * diff(B.zS_bnd);
        S.svp  = B.svp_bnd(1)  + lhs(:,3) * diff(B.svp_bnd);

    case 'normal'
        S.freq = freq0 + B.sig_freq * norminv(lhs(:,1));
        S.zS   = zS0   + B.sig_zS   * norminv(lhs(:,2));
        S.svp  =         B.sig_svp  * norminv(lhs(:,3));   % centered at 0

    otherwise
        error('Unknown distribution type: %s', dist.type);
end
end
