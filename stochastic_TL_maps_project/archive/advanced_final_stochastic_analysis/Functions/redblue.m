function cmap = redblue(n)
% Red-white-blue diverging colormap (n steps, default 256).
    if nargin < 1, n = 256; end
    h = n / 2;
    r1 = [linspace(0.2,1,h); linspace(0.2,1,h); ones(1,h)]';
    r2 = [ones(1,h); linspace(1,0.2,h); linspace(1,0.2,h)]';
    cmap = [r1; r2];
end
