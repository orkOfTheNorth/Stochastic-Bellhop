function p = jensenHalfspaceImage(freq, zS, r, z, c)
% jensenHalfspaceImage  Point source in a fluid half-space — METHOD OF
% IMAGES ONLY. No normal modes, no waveguide, no bottom reflection.
%
% Source: Jensen, Kuperman, Porter, Schmidt (2011), "Computational Ocean
% Acoustics", 2nd ed., Sec. 2.3.4 "Point Source in Fluid Halfspace", p.80,
% Eq. (2.76):
%
%   psi(r,z) = S_w * [ exp(ik R1)/R1  -  exp(ik R2)/R2 ]
%
%   R1 = sqrt(r^2 + (z - zS)^2)   real source at depth zS
%   R2 = sqrt(r^2 + (z + zS)^2)   image source at depth -zS
%
% The minus sign between the two terms is exactly the pressure-release
% (free) surface boundary condition at z=0 (Jensen eq. 2.69/2.73): the
% image source has opposite sign so p(r,0) = 0 identically for any r.
%
% The 1/(4*pi) normalization constant in Jensen's eq. (2.73)/(2.76) is
% dropped here (S_w = 1, no 4*pi) because it cancels exactly in the TL
% ratio TL = -20*log10(|p|/|p_ref|) as long as p_ref uses the same
% convention (see pressureToTL.m and its use in
% run_jensen_halfspace_validation.m, where p_ref = 1 at R=1m).
%
% Inputs:
%   freq — source frequency [Hz]
%   zS   — source depth [m]
%   r, z — range [m] and depth [m] grids (same size), e.g. from meshgrid
%   c    — sound speed [m/s] (isovelocity — this is a HOMOGENEOUS half-space)
%
% Output:
%   p — complex pressure amplitude at each (r,z), same size as r/z

k = 2*pi*freq / c;

R1 = sqrt(r.^2 + (z - zS).^2);
R2 = sqrt(r.^2 + (z + zS).^2);

p = exp(1i*k*R1)./R1 - exp(1i*k*R2)./R2;
end
