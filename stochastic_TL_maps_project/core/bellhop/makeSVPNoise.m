function svp = makeSVPNoise(svp_level, max_depth, cfg)
% Generate a perturbed Sound Velocity Profile (SVP) using coupled T_surface + MLD.
%
% The single parameter svp_level (degC) simultaneously shifts:
%   T_surface = T_NOM + svp_level
%   MLD       = MLD_NOM + MLD_COEFF * svp_level
%
% Physical validity:
%   T_surface clamped to [10, 35] degC
%   MLD       clamped to [5 m, 0.8 * max_depth]
%
% Inputs:
%   svp_level — signed scalar (degC).  Positive = warmer surface + deeper MLD.
%   max_depth — water column depth (m)
%   cfg       — (optional) struct from loadConfig(); uses defaults if omitted.
%
% Output:
%   svp — (200 x 2) matrix: [depth_m, sound_speed_m_s]

if nargin < 3 || isempty(cfg)
    T_NOM     = 25;
    MLD_NOM   = 30;
    MLD_COEFF = 5;
    S_sal     = 37;
else
    T_NOM     = cfg.SVP.nominal_surface_temp_C;
    MLD_NOM   = cfg.SVP.nominal_MLD_m;
    MLD_COEFF = cfg.SVP.MLD_coupling_m_per_C;
    S_sal     = cfg.SVP.salinity_ppt;
end

T_surf = max(10, min(35,             T_NOM + svp_level));
MLD    = max(5,  min(0.8*max_depth,  MLD_NOM + MLD_COEFF * svp_level));

depths = linspace(0, max_depth, 200);

d_thermo_bot = min(180, max_depth);
d_deep       = min(400, max_depth);
ctrl_d = unique([0; MLD; d_thermo_bot; d_deep; max_depth]);

base_d = [0;  30;  180; 400;  5000];
base_t = [T_surf; T_surf; 17; 13.6; 13.6];
ctrl_t = interp1(base_d, base_t, ctrl_d, 'linear', 'extrap');
ctrl_t(ctrl_d <= MLD) = T_surf;   % isothermal mixed layer

t  = interp1(ctrl_d, ctrl_t, depths, 'linear', 'extrap');
sv = 1499.2 + 4.6*t - 0.055*t.^2 + 0.00029*t.^3 + ...
     (1.34 - 0.01*t).*(S_sal - 35) + 0.016*depths;

svp = [depths.' sv.'];
end
