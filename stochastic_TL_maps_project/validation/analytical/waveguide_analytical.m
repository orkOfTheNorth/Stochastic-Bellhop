function TL = waveguide_analytical(r_vec, z_vec, freq, c, D, zS, alpha_dBpkm)
% waveguide_analytical  Isovelocity waveguide TL via normal modes.
%
% Rigid bottom (Neumann), pressure-release surface (Dirichlet).
% Includes 1/sqrt(r) cylindrical spreading (matches Bellhop/KRAKEN convention).
%
% Mode shapes: k_zm = (m+0.5)*pi/D
% Pressure:    p(r,z) = sqrt(pi/2) * sum_m  Nm * sin(k_zm*zs)*sin(k_zm*z)
%                        * exp(-1i*(k_xm*r - pi/4)) / (D * sqrt(k_xm*r))
% TL referenced to point source at 1 m.
%
% Inputs:
%   r_vec       — range vector [1 × Nr], metres  (avoid r=0)
%   z_vec       — depth vector, metres
%   freq        — frequency (Hz)
%   c           — sound speed (m/s)
%   D           — water depth (m)
%   zS          — source depth (m)
%   alpha_dBpkm — Thorp absorption (dB/km); use 0 for lossless

r_vec = r_vec(:)';          % [1 × Nr]
z_vec = z_vec(:);           % [Nz × 1]

k     = 2 * pi * freq / c;
% Convert absorption: dB/km → nepers/m
alpha_npm = alpha_dBpkm / (20 * log10(exp(1)) * 1000);

P = zeros(numel(z_vec), numel(r_vec));

for m = 0:3000
    k_zm   = (m + 0.5) * pi / D;
    inside = k^2 - k_zm^2;
    if inside <= 0
        break;   % all remaining modes evanescent
    end
    k_xm   = sqrt(inside) - 1i * alpha_npm;

    mode_zS = sin(k_zm * zS);
    if abs(mode_zS) < 1e-12, continue; end

    mode_z = sin(k_zm .* z_vec);   % [Nz × 1]

    % Cylindrical spreading: 1/sqrt(k_xm * r), phase correction -pi/4
    prop = exp(-1i * (k_xm .* r_vec - pi/4)) ./ sqrt(real(k_xm) .* r_vec);

    % Normalization following KRAKEN convention
    Nm   = sqrt(pi / 2) / D;
    term = Nm * mode_zS .* (mode_z * prop);
    P    = P + term;
end

TL = -20 * log10(max(abs(P), 1e-20));
end
