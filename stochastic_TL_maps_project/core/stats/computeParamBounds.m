function B = computeVarianceBounds(cfg, dist)
% Compute perturbation bounds and variances for one distribution config entry.
%
% Inputs:
%   cfg  — struct from loadConfig()
%   dist — one element of cfg.distributions (has .type and .level_pct)
%
% Output B (struct):
%   B.freq_bnd   [lo, hi] Hz
%   B.zS_bnd     [lo, hi] m
%   B.svp_bnd    [lo, hi] deg-C  (in svp_level units)
%   B.var_freq   Hz^2
%   B.var_zS     m^2
%   B.var_svp    degC^2
%   B.sig_freq   Hz   (std dev, for LHS Normal sampling)
%   B.sig_zS     m
%   B.sig_svp    degC

freq0 = cfg.nominal.freq_Hz;
zS0   = cfg.nominal.zS_m;
MLD0  = cfg.SVP.nominal_MLD_m;
kMLD  = cfg.SVP.MLD_coupling_m_per_C;   % m / degC
x     = dist.level_pct;                 % error level (%)

switch lower(dist.type)
    case 'uniform'
        % Bound = nominal * x/100; 2*bound is total range
        bnd_freq = freq0 * x / 100;
        bnd_zS   = zS0   * x / 100;
        bnd_MLD  = MLD0  * x / 100;          % MLD shift bound in metres
        bnd_svp  = bnd_MLD / kMLD;           % convert to degC

        B.freq_bnd = [freq0 - bnd_freq, freq0 + bnd_freq];
        B.zS_bnd   = [zS0   - bnd_zS,   zS0   + bnd_zS];
        B.svp_bnd  = [-bnd_svp,          bnd_svp];

        B.var_freq = (2*bnd_freq)^2 / 12;
        B.var_zS   = (2*bnd_zS)^2   / 12;
        B.var_svp  = (2*bnd_svp)^2  / 12;

        B.sig_freq = sqrt(B.var_freq);
        B.sig_zS   = sqrt(B.var_zS);
        B.sig_svp  = sqrt(B.var_svp);

    case 'normal'
        % 2*sigma = x% of nominal
        B.sig_freq = freq0 * x / 200;
        B.sig_zS   = zS0   * x / 200;
        sig_MLD    = MLD0  * x / 200;
        B.sig_svp  = sig_MLD / kMLD;

        % Bounds = ±3*sigma for practical sampling range
        B.freq_bnd = [freq0 - 3*B.sig_freq, freq0 + 3*B.sig_freq];
        B.zS_bnd   = [zS0   - 3*B.sig_zS,   zS0   + 3*B.sig_zS];
        B.svp_bnd  = [-3*B.sig_svp,           3*B.sig_svp];

        B.var_freq = B.sig_freq^2;
        B.var_zS   = B.sig_zS^2;
        B.var_svp  = B.sig_svp^2;

    otherwise
        error('Unknown distribution type: %s', dist.type);
end
end
