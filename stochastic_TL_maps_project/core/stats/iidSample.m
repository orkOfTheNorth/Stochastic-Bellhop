function S = iidSample(N, cfg, dist, B, seed)
% Generate N IID (pure random) samples — alternative to LHS stratification.
% Same interface as lhsSample but uses randn/rand instead of lhsdesign.
if nargin < 5, seed = cfg.MC.rng_seed + 1000; end
rng(seed);
freq0 = cfg.nominal.freq_Hz;
zS0   = cfg.nominal.zS_m;
switch lower(dist.type)
    case 'normal'
        S.freq = freq0 + B.sig_freq * randn(N, 1);
        S.zS   = zS0   + B.sig_zS   * randn(N, 1);
        S.svp  =         B.sig_svp  * randn(N, 1);
    case 'uniform'
        S.freq = B.freq_bnd(1) + rand(N,1) * diff(B.freq_bnd);
        S.zS   = B.zS_bnd(1)   + rand(N,1) * diff(B.zS_bnd);
        S.svp  = B.svp_bnd(1)  + rand(N,1) * diff(B.svp_bnd);
    otherwise
        error('Unknown distribution type: %s', dist.type);
end
end
