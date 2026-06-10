function out = lpf2d(img, cutoff_z, cutoff_r)
% Ideal rectangular 2-D low-pass filter in the FFT domain.
%
%   cutoff_z — fraction of z-axis (depth) bandwidth to KEEP  [0,1]
%   cutoff_r — fraction of r-axis (range) bandwidth to KEEP  [0,1]
%              If omitted, cutoff_r = cutoff_z  (isotropic)
%
%   Examples (keeping = blocks features SMALLER than):
%     cutoff = 0.10 → keeps only the lowest 10% of frequencies (strong filter)
%     cutoff = 0.70 → keeps lowest 70% (mild filter)
%
%   Physical helper:
%     cutoff = 2 * pixel_spacing / min_feature_size
%
% Output is clamped to [0, 1] (suitable for probability maps).

if nargin < 2 || isempty(cutoff_z)
    cutoff_z = 0.90;
end
if nargin < 3 || isempty(cutoff_r)
    cutoff_r = cutoff_z;
end

[Nz, Nr] = size(img);
F = fftshift(fft2(double(img)));

fz = (-Nz/2 : Nz/2-1) / Nz;
fr = (-Nr/2 : Nr/2-1) / Nr;
[FR, FZ] = meshgrid(fr, fz);

mask = (abs(FZ) <= cutoff_z/2) & (abs(FR) <= cutoff_r/2);
out  = real(ifft2(ifftshift(F .* mask)));
out  = max(0, min(1, out));
end
