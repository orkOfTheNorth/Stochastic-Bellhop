function TL = pressureToTL(p, p_ref)
% pressureToTL  Super-basic pressure-to-TL conversion.
%
%   TL = -20*log10(|p| / |p_ref|)      [dB]
%
% p_ref defaults to 1 (a source normalized to unit pressure amplitude at
% 1 m in free space — Jensen (2011) Computational Ocean Acoustics, eq. 2.77-2.78).
%
% Inputs:
%   p     — complex (or real) pressure amplitude, any size
%   p_ref — reference pressure (default 1)
%
% Output:
%   TL — transmission loss in dB, same size as p

if nargin < 2 || isempty(p_ref), p_ref = 1; end
TL = -20 * log10(abs(p) ./ abs(p_ref));
end
