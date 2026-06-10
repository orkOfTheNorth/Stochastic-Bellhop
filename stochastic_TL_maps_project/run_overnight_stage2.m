try
  cd('c:\Users\orind\Stochastic-Bellhop\stochastic_TL_maps_project');
  addpath(genpath('Shared_Utils'));
  addpath(genpath('Bellhop'));
  disp('=== Stage 2: run_Comparison ===');
  run_Comparison;
  disp('=== Stage 2: run_LN3 ===');
  run_LN3;
  disp('=== Stage 2: re-running run_pce (for new combos) ===');
  run_pce;
  disp('=== Stage 2 complete ===');
catch ME
  fid = fopen('overnight_stage2_error.log','w');
  fprintf(fid, '%s\n', ME.message);
  fclose(fid);
end
exit;
