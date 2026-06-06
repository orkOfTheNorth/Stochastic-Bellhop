function progressBar(i, N, label, t_start)
% Print a visual ASCII progress bar to the command window.
%
% i       : current iteration (1-based)
% N       : total iterations
% label   : short string shown before the bar (e.g. subset name)
% t_start : tic value from start of loop (for ETA)
%
% Example output:
%   [zS_freq_svp]  [████████████░░░░░░░░░░░░]  24/50 (48%)  elapsed 72s  ETA 78s

    BAR_LEN = 24;
    filled  = round(BAR_LEN * i / N);
    empty   = BAR_LEN - filled;
    bar     = [repmat(char(9608), 1, filled), repmat(char(9617), 1, empty)];

    elapsed = toc(t_start);
    if i > 0
        eta = elapsed / i * (N - i);
    else
        eta = 0;
    end

    fprintf('  [%-14s]  [%s]  %d/%d (%3.0f%%)  elapsed %ds  ETA %ds\n', ...
            label, bar, i, N, 100*i/N, round(elapsed), round(eta));
end
