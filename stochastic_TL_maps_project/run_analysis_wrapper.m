try
    run_analysis_only
catch ME
    fprintf('ERROR: %s\n', ME.message);
    for k = 1:numel(ME.stack)
        fprintf('  at %s line %d\n', ME.stack(k).file, ME.stack(k).line);
    end
    exit(1)
end
