function saveSVPGraph(svp, save_path, label)
% Save a sound-speed vs depth PNG for visual validation of SVP perturbation.
%
% Inputs:
%   svp       — Nx2 matrix [depth_m, speed_m_s]
%   save_path — full path WITHOUT extension
%   label     — title string (e.g. 'svp_level=+0.30 degC')

fig = figure('Visible', 'off', 'Position', [50 50 420 580]);
plot(svp(:,2), svp(:,1), 'b-', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse');
xlabel('Sound Speed (m/s)');
ylabel('Depth (m)');
title(label, 'Interpreter', 'none');
grid on;
saveFigPNG(fig, save_path);
close(fig);
end
