%% השוואת אי-ודאות: מונטה קרלו מול שיטת הדלתא
clear; close all; clc;

% פרמטרים קבועים
f = 10000;          % 10 kHz
rmax = 50000;       % 50 km
z_src = 4.5;         % עומק מקור (הנחה)
svp = 'summer';     % פרופיל קיץ
fom = 0;            % FOM 0dB

% הגדרת אי-ודאות עבור עומק הקרקעית (Uniform [5, 100])
z_min = 5;
z_max = 100;
z_nominal = (z_min + z_max) / 2; % 52.5 m
sigma_input_sq = (z_max - z_min)^2 / 12; % שונות של התפלגות אחידה

%% 1. שיטת ה-Delta Method (חישוב מקומי)
fprintf('Calculating Delta Method (Local Gradient)...\n');
eps = 0.5; % "תזוזה" קטנה לחישוב נגזרת

% הרצת בסיס והרצת נגזרת
sim_base = {f, rmax, z_src, 0, 0, svp, sprintf('const_%d', round(z_nominal)), [1600, 1.8, 0.5], fom};
sim_nudge = {f, rmax, z_src, 0, 0, svp, sprintf('const_%d', round(z_nominal + eps)), [1600, 1.8, 0.5], fom};

TL_base = simpleBellhopHazat(sim_base);
TL_nudge = simpleBellhopHazat(sim_nudge);

% חישוב שונות לפי Delta Method
gradient = (TL_nudge - TL_base) / eps;
Var_Delta = (gradient.^2) * sigma_input_sq;

%% 2. שיטת Monte Carlo (דגימה גלובלית)
N_runs = 30; % מספר ריצות
fprintf('Running Monte Carlo (%d iterations)...\n', N_runs);

[rows, cols] = size(TL_base);
TL_MC_stack = zeros(rows, cols, N_runs);
z_samples = z_min + (z_max - z_min) * rand(N_runs, 1); % הגרלה יוניפורמית

for i = 1:N_runs
    sim_mc = {f, rmax, z_src, 0, 0, svp, sprintf('const_%d', round(z_samples(i))), [1600, 1.8, 0.5], fom};
    TL_MC_stack(:,:,i) = simpleBellhopHazat(sim_mc);
    if mod(i,10)==0, fprintf('Iteration %d/%d\n', i, N_runs); end
end

Var_MC = var(TL_MC_stack, 0, 3);

%% 3. ויזואליזציה והשוואה
ranges = linspace(0, rmax, cols);
depths = linspace(0, z_max, rows); % תצוגה עד עומק הקרקעית המקסימלי

fig = figure('Color', 'w', 'Position', [100 100 1200 500]);

% מפת שונות - מונטה קרלו
ax1 = subplot(1,2,1);
imagesc(ranges/1000, depths, Var_MC);
title('Monte Carlo Variance (Global)');
xlabel('Range (km)'); ylabel('Depth (m)');
colorbar; colormap(ax1, 'hot'); clim([0 max(Var_MC(:))*0.5]);

% מפת שונות - שיטת הדלתא
ax2 = subplot(1,2,2);
imagesc(ranges/1000, depths, Var_Delta);
title('Delta Method Variance (Local)');
xlabel('Range (km)'); ylabel('Depth (m)');
colorbar; colormap(ax2, 'hot'); clim([0 max(Var_MC(:))*0.5]);

% הוספת Tooltip אינטראקטיבי להצגת השונות בנקודה
dcm = datacursormode(fig);
dcm.Enable = 'on';
dcm.UpdateFcn = @(obj, event_obj) sprintf('Range: %.1f km\nDepth: %.1f m\nVar MC: %.2f\nVar Delta: %.2f', ...
    event_obj.Position(1), event_obj.Position(2), ...
    interp2(ranges/1000, depths, Var_MC, event_obj.Position(1), event_obj.Position(2)), ...
    interp2(ranges/1000, depths, Var_Delta, event_obj.Position(1), event_obj.Position(2)));

linkaxes([ax1, ax2], 'xy');