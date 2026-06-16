function progressBar(i, N, label, t_start)
% ASCII progress bar with ETA printed to the command window.
%
% Inputs:
%   i       — current iteration (1-based)
%   N       — total iterations
%   label   — short label string (truncated to 14 chars for alignment)
%   t_start — tic value from loop start

BAR_LEN = 24;
frac    = i / N;
filled  = round(frac * BAR_LEN);
bar     = [repmat(char(9608), 1, filled), repmat(char(9617), 1, BAR_LEN - filled)];

elapsed = toc(t_start);
if i > 0
    eta = elapsed / i * (N - i);
else
    eta = 0;
end

lbl = sprintf('%-14s', label);
fprintf('[%s]  [%s]  %d/%d (%.0f%%)  elapsed %.0fs  ETA %.0fs\n', ...
        lbl, bar, i, N, frac*100, elapsed, eta);
end
