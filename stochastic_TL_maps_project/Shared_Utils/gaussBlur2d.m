function out = gaussBlur2d(img, sigma_z_px, sigma_r_px)
% Separable 2-D Gaussian blur with potentially different sigma per axis.
%
% Inputs:
%   img        — 2-D matrix to blur
%   sigma_z_px — std dev in DEPTH (row) direction, in pixels
%   sigma_r_px — std dev in RANGE (col) direction, in pixels
%
% Skips convolution on an axis when sigma < 0.5 px (negligible effect).
% Output is clamped to [0, 1].

img = double(img);

if sigma_z_px >= 0.5
    rz  = ceil(3 * sigma_z_px);
    kz  = exp(-((-rz:rz).^2) / (2 * sigma_z_px^2));
    kz  = kz / sum(kz);
    img = conv2(kz(:), 1, img, 'same');
end

if sigma_r_px >= 0.5
    rr  = ceil(3 * sigma_r_px);
    kr  = exp(-((-rr:rr).^2) / (2 * sigma_r_px^2));
    kr  = kr / sum(kr);
    img = conv2(1, kr(:)', img, 'same');
end

out = max(0, min(1, img));
end
