function [names, labels, active] = subsetDefs()
% Canonical 7 uncertainty subsets — single source of truth.
%
% active{s} = [a_freq, a_zS, a_svp]  (1 = parameter included, 0 = fixed at nominal)
%
% Subset order:
%   1 zS          [0 1 0]
%   2 freq        [1 0 0]
%   3 svp         [0 0 1]
%   4 zS_freq     [1 1 0]
%   5 zS_svp      [0 1 1]
%   6 freq_svp    [1 0 1]
%   7 zS_freq_svp [1 1 1]

names = {'zS','freq','svp','zS_freq','zS_svp','freq_svp','zS_freq_svp'};

labels = {'z_S only', 'Freq only', 'SVP only', ...
          'z_S+Freq', 'z_S+SVP',  'Freq+SVP', 'z_S+Freq+SVP'};

active = {[0 1 0], [1 0 0], [0 0 1], ...
          [1 1 0], [0 1 1], [1 0 1], [1 1 1]};
end
