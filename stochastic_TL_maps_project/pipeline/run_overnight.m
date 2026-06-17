try
  cd('c:\Users\orind\Stochastic-Bellhop\stochastic_TL_maps_project');
  addpath(genpath(fullfile(ROOT, 'core')));
addpath(fullfile(ROOT, 'binaries'));
    run_MC;
  run_pce;
catch ME
  fid = fopen('overnight_error.log','w');
  fprintf(fid, '%s\n', ME.message);
  fclose(fid);
end
exit;
