function done = skipIfDone(sentinel, label)
% Returns true and prints a skip message if the sentinel file exists.
% Use at the top of each pipeline stage to avoid re-running completed work.
%
%   if skipIfDone(fullfile(out_dir,'result.mat'), 'baseline/Normal_5pct')
%       continue;
%   end

done = isfile(sentinel);
if done
    fprintf('  SKIP (done): %s\n', label);
end
end
