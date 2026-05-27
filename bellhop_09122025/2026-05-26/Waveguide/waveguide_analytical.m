function TL = waveguide_analytical(r_vec, z_vec, freq, c, D, zS, alpha_dBpkm)
%% waveguide_analytical  –  Ideal waveguide TL via normal mode sum
%
% Boundary conditions (matching the simple analytical solution):
%   Top  : pressure-release   (P = 0 at z = 0)
%   Bottom: rigid             (dP/dz = 0 at z = D)
%
% Mode eigenvalues:  k_zm = (m + 0.5) * pi / D,  m = 0, 1, 2, ...
%
% This is the 3D POINT SOURCE version (cylindrical spreading).
% The far-field Hankel function approximation is used:
%   H0(1)(k_xm * r) ≈ sqrt(2 / (pi * k_xm * r)) * exp(j*(k_xm*r - pi/4))
% This adds the 1/sqrt(r) cylindrical spreading that Bellhop also uses.
%
% Inputs
%   r_vec        range vector  [m]
%   z_vec        depth vector  [m]
%   freq         frequency     [Hz]
%   c            sound speed   [m/s]  constant SVP
%   D            water depth   [m]
%   zS           source depth  [m]
%   alpha_dBpkm  volume attenuation [dB/km], 0 = no absorption
%
% Output
%   TL   transmission loss  [dB],  size Nz × Nr

% ── wavenumber (complex if volume absorption > 0) ────────────────────────
k = 2*pi*freq/c;
alpha_npm = alpha_dBpkm * log(10) / (20000);   % dB/km → Np/m
k = k + 1i*alpha_npm;

r_vec = r_vec(:)';     % 1 × Nr
z_vec = z_vec(:);      % Nz × 1

% Avoid division by zero at r=0
r_vec(r_vec < 1e-3) = 1e-3;

P = zeros(numel(z_vec), numel(r_vec));

% ── mode sum ─────────────────────────────────────────────────────────────
m = 0;
while true
    k_zm = (m + 0.5) * pi / D;

    % Stop once we reach the evanescent regime.
    % Evanescent modes die out within < 1 m at 10 kHz — safe to ignore.
    if k_zm > real(k)
        break
    end

    k_xm = sqrt(k^2 - k_zm^2);   % horizontal wavenumber (complex ok)

    % Far-field Hankel function: cylindrical spreading + phase
    % H0(1)(k_xm*r) ≈ sqrt(2/(pi*k_xm*r)) * exp(j*(k_xm*r - pi/4))
    mode_x = sqrt(2 ./ (pi * k_xm * r_vec)) .* exp(1i*(k_xm*r_vec - pi/4));

    % Mode shapes at source and receiver
    mode_src = sin(k_zm * zS);     % scalar
    mode_rec = sin(k_zm * z_vec);  % Nz × 1

    % Amplitude factor (from the modal Green's function formula)
    amp = 1 / (k_xm * D);

    % Add mode contribution: outer product gives (Nz × Nr)
    P = P + amp * mode_src * (mode_rec * mode_x);

    m = m + 1;
end

% ── transmission loss ─────────────────────────────────────────────────────
% The absolute level of this formula has an unknown normalisation constant
% relative to Bellhop's reference (1 Pa at 1 m free-field).
% The constant is factored out in run_waveguide_comparison.m by calibrating
% both models to the same depth-averaged level at a reference range.
TL = -20*log10(abs(P) + eps);

end
