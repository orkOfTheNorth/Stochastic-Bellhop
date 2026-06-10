function cmap = redblue(n)
% Red-white-blue diverging colormap.
%
%   Lower half: blue  [0.2, 0.2, 1.0] → white [1, 1, 1]
%   Upper half: white [1, 1, 1]       → red   [1, 0.2, 0.2]

if nargin < 1 || isempty(n)
    n = 256;
end

half = floor(n / 2);
lo   = [linspace(0.2, 1, half)', linspace(0.2, 1, half)', ones(half, 1)];
hi   = [ones(n-half, 1), linspace(1, 0.2, n-half)', linspace(1, 0.2, n-half)'];
cmap = [lo; hi];
end
