function svp = makeSVPNoise(svp_level, max_depth)
% Unified SVP noise parameter: simultaneously shifts surface temperature
% AND mixing layer depth (MLD) using an empirical physical coupling.
%
% svp_level : signed scalar in °C (positive = warmer surface + deeper MLD)
% max_depth : water column depth in metres (SVP is built to this depth)
%
% Physical model (Mediterranean summer):
%   T_surface = 25 + svp_level        [nominal 25°C]
%   MLD       = 30 + 5*svp_level      [nominal 30m, 5 m/°C coupling]
%
% Physical validity limits:
%   T_surface clamped to [10, 35]°C
%   MLD       clamped to [5m, 0.8*max_depth]
%
% For distribution study error levels:
%   1%  → svp_level ∈ [−0.25, +0.25]°C  → ΔT=0.25°C,  ΔMLD=1.25m
%   5%  → svp_level ∈ [−1.25, +1.25]°C  → ΔT=1.25°C,  ΔMLD=6.25m
%   10% → svp_level ∈ [−2.50, +2.50]°C  → ΔT=2.50°C,  ΔMLD=12.5m

    if nargin < 2, max_depth = 35; end

    MLD_NOM   = 30;   % m
    T_NOM     = 25;   % °C
    MLD_COEFF = 5;    % m / °C

    T_surf = max(10, min(35,   T_NOM + svp_level));
    MLD    = max(5,  min(0.8 * max_depth, MLD_NOM + MLD_COEFF * svp_level));

    depths = linspace(0, max_depth, 200);

    % Build temperature profile control points
    d_thermo_bot = min(180, max_depth);
    d_deep       = min(400, max_depth);
    ctrl_d = unique([0; MLD; d_thermo_bot; d_deep; max_depth]);
    % Interpolate from base summer profile to get temps at those control points
    base_d = [0;  30;  180; 400;  5000];
    base_t = [T_surf; T_surf; 17; 13.6; 13.6];  % surface already shifted
    ctrl_t = interp1(base_d, base_t, ctrl_d, 'linear', 'extrap');
    % Override isothermal layer up to MLD
    ctrl_t(ctrl_d <= MLD) = T_surf;

    t  = interp1(ctrl_d, ctrl_t, depths, 'linear', 'extrap');
    S  = 37;   % salinity [ppt]
    sv = 1499.2 + 4.6*t - 0.055*t.^2 + 0.00029*t.^3 + ...
         (1.34 - 0.01*t).*(S - 35) + 0.016*depths;
    svp = [depths.' sv.'];
end
