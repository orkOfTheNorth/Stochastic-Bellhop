function createBellhopBty(bathymetry) % Create bellhop.bty
    fileID = fopen('bellhop.bty', 'w');
    fprintf(fileID, ["'L'"]);
    fprintf(fileID, ['\n']);
    fprintf(fileID, [num2str(size(bathymetry, 1)) '\n']);
    for i = 1:size(bathymetry,1)
      r = bathymetry(i,1)/1000;
      d = bathymetry(i,2);
      fprintf(fileID, [num2str(r) '   ' num2str(d) '   /\n']);
    end
    fclose(fileID);
end