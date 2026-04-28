%% 2) Monte Carlo Method for TL Variance (Fixed Grid & NaN Handling)
clear; close all; clc; warning('off');
addpath('../Code/Functions'); % נתיב מעודכן

%% User Parameters & Uncertainty Definitions
sourceFrequency = 10000;    % Hz
maxRange = 50000;           % m
sourceDepth = 4.5;          % m
sourceHalfBeam = 0;         % deg
sourceTilt = 0;             % deg
svpType = "summer";         % string
geoAcoustics = [1500 1.63 0.07]; % [c_bottom rho alpha]
fom = 0;                    % dB

% Monte Carlo Parameters
N_runs = 50;
z_min = 5;
z_max = 100;

%% 1. Define Master Grid
% רשת אחידה שאליה כל הריצות יעברו אינטרפולציה
z_master = linspace(0, z_max, 500); % עומק מ-0 עד 100 קבוע
r_master = linspace(0, maxRange, 1000); % טווח קבוע
[R_mast, Z_mast] = meshgrid(r_master, z_master);

%% 2. Generate Spaced-Out Samples (Latin Hypercube)
rng('default'); 
normalized_samples = lhsdesign(N_runs, 1); 
z_samples = z_min + (z_max - z_min) * normalized_samples; 
z_samples = z_samples(randperm(N_runs));

%% 3. Initialize Storage
% עכשיו נאתחל את המטריצה התלת-מימדית מראש לפי גודל ה-Master Grid
TL_all = NaN(length(z_master), length(r_master), N_runs);

%% 4. Execute Monte Carlo Loop
for i = 1:N_runs
    fprintf('Run %d/%d: z = %.2f m...\n', i, N_runs, z_samples(i));
    bathy_str = sprintf('const_%.2f', z_samples(i));
    sim_pars = {sourceFrequency, maxRange, sourceDepth, sourceHalfBeam, ...
        sourceTilt, svpType, bathy_str, geoAcoustics, fom};
    
    % הרצת המודל
    [TL_current, r_current, z_current] = simpleBellhopHazat(sim_pars);
    
    % ניקוי נתונים: הפיכת ערכי TL קיצוניים או שגיאות ל-NaN
    % (ערכים מעל 150dB הם לרוב רעש אקוסטי אפסי / "מוות" של קרן)
    TL_current(TL_current > 150) = NaN;
    TL_current(isinf(TL_current)) = NaN;
    TL_current(TL_current == 0) = NaN;
    
    % יצירת רשת זמנית עבור הריצה הנוכחית
    [R_curr, Z_curr] = meshgrid(r_current, z_current);
    
    % אינטרפולציה לרשת הקבועה (ערכים מחוץ לתחום יהפכו ל-NaN)
    TL_interp = interp2(R_curr, Z_curr, TL_current, R_mast, Z_mast, 'linear', NaN);
    
    % שמירה למערך התלת-מימדי
    TL_all(:,:,i) = TL_interp;
end

%% 5. Calculate Statistics (Ignoring NaNs)
% שימוש ב-'omitnan' כדי שהשונות לא תושפע מאזורים "מתים" בריצות רדודות
E_TL_mc = mean(TL_all, 3, 'omitnan');
Var_TL_mc = var(TL_all, 0, 3, 'omitnan');

% חזרה לשמות המשתנים המקוריים בשביל המשך התאימות לקודים האחרים
r_grid = r_master;
z_grid = z_master;

%% 6. Save Data for Comparison Script
save('monte_carlo_results.mat', 'r_grid', 'z_grid', 'E_TL_mc', 'Var_TL_mc', 'N_runs', 'z_samples');
fprintf('Results saved to monte_carlo_results.mat\n');

%% Visualization with Interactive Hover
fig = figure('Name', 'Monte Carlo: Expected Value and Variance', 'Position', [150, 150, 1200, 500]);

% Plot Expected Value
subplot(1,2,1);
imagesc(r_grid/1000, z_grid, E_TL_mc);
colorbar; set(gca, 'YDir', 'reverse');
title(sprintf('Expected Value E[TL] (%d MC Runs)', N_runs));
xlabel('Range (km)'); ylabel('Depth (m)');
colormap('jet');

% Plot Variance
subplot(1,2,2);
imagesc(r_grid/1000, z_grid, Var_TL_mc);
colorbar; set(gca, 'YDir', 'reverse');
title('Variance Var(TL) via Monte Carlo');
xlabel('Range (km)'); ylabel('Depth (m)');
colormap('hot');
% נגביל את התצוגה של השונות כדי שרעשי קצה קטנים לא ישטפו את הצבע
clim_max = prctile(Var_TL_mc(:), 95); % חיתוך צבע אחוזון 95
if clim_max > 0 && ~isnan(clim_max)
    clim([0 clim_max]);
end

% Set up custom data cursor
dcm_obj = datacursormode(fig);
set(dcm_obj, 'UpdateFcn', {@customDatatip, r_grid, z_grid, E_TL_mc, Var_TL_mc});

%% Custom Datatip Function
function txt = customDatatip(~, event_obj, r_grid, z_grid, E_TL, Var_TL)
    pos = get(event_obj, 'Position');
    clicked_r = pos(1) * 1000; 
    clicked_z = pos(2);
    
    [~, r_idx] = min(abs(r_grid - clicked_r));
    [~, z_idx] = min(abs(z_grid - clicked_z));
    
    txt = {sprintf('Range: %.2f km', pos(1)), ...
           sprintf('Depth: %.1f m', pos(2)), ...
           sprintf('E[TL]: %.2f dB', E_TL(z_idx, r_idx)), ...
           sprintf('Var(TL): %.2f dB^2', Var_TL(z_idx, r_idx))};
end