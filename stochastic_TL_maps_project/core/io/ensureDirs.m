function ensureDirs(varargin)
% ensureDirs(dir1, dir2, ...) — create each directory path if it does not exist.
for k = 1:nargin
    d = varargin{k};
    if ~exist(d, 'dir')
        mkdir(d);
    end
end
end
