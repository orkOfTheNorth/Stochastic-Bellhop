function saveFigPNG(fig, base_path)
% Save figure as a 300-dpi PNG. No PDF output.
%
% Inputs:
%   fig       — figure handle
%   base_path — full path WITHOUT extension (extension .png is added)

png_path = [base_path '.png'];
exportgraphics(fig, png_path, 'Resolution', 300);
end
