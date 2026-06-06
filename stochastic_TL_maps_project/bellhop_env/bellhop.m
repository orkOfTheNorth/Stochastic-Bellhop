function bellhop( filename )
%filename
%typeinfo(filename)
% run the BELLHOP program
%
% usage: bellhop( filename )
% where filename is the environmental file

runbellhop = which( 'bellhop.exe' );
% run
%typeinfo(runbellhop)
if ( isempty( runbellhop ) )
   error( 'bellhop.exe not found in your Octave path' )
else
   system( [ '"' runbellhop '" ' filename ] );
end