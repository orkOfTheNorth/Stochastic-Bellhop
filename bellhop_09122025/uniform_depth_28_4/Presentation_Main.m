%% Presentation Main - Bellhop Uncertainty Analysis
% סקריפט לריכוז והצגת כלל התוצאות לפרזנטציה
clear; close all; clc; warning('off');
addpath('../Code/Functions'); % קריטי להרצות הבסיס שנוספו

% בדיקת קיום קבצים
files = {'monte_carlo_results.mat', 'delta_method_results.mat'};
for f = 1:length(files)
    if ~exist(files{f}, 'file')
        error('חסר קובץ נתונים: %s. וודא שהרצת את כל הסימולציות.', files{f});
    end
end

%% 0. הרצות בסיס חיות (Baseline Runs: 52.5m & 100m)
fprintf('מריץ סימולציות בסיס (52.5m, 100m) עבור הפרזנטציה...\n');
sourceFreq = 10000; maxR = 50000; srcDepth = 4.5;
geoAc = [1500 1.63 0.07]; fom = 0; svp = "summer";

% הרצה 52.5 מטר
sim_pars_52 = {sourceFreq, maxR, srcDepth, 0, 0, svp, 'const_52.5', geoAc, fom};
[TL_52, r_52, z_52] = simpleBellhopHazat(sim_pars_52);
TL_52(TL_52 > 150 | TL_52 == 0 | isinf(TL_52)) = NaN; % ניקוי רעשים לתצוגה

% הרצה 100 מטר
sim_pars_100 = {sourceFreq, maxR, srcDepth, 0, 0, svp, 'const_100.0', geoAc, fom};
[TL_100, r_100, z_100] = simpleBellhopHazat(sim_pars_100);
TL_100(TL_100 > 150 | TL_100 == 0 | isinf(TL_100)) = NaN; % ניקוי רעשים לתצוגה

fig0 = figure('Name', '0. Baseline Runs', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.3 0.3 0.4 0.4]);
tiledlayout(1,2);

nexttile;
imagesc(r_52/1000, z_52, TL_52);
title('Baseline TL: z = 52.5 m (Nominal)'); colorbar; colormap(gca, 'jet');
set(gca, 'YDir', 'reverse'); xlabel('Range (km)'); ylabel('Depth (m)');

nexttile;
imagesc(r_100/1000, z_100, TL_100);
title('Baseline TL: z = 100.0 m (Max Depth)'); colorbar; colormap(gca, 'jet');
set(gca, 'YDir', 'reverse'); xlabel('Range (km)'); ylabel('Depth (m)');

%% 1. תצוגת מונטה קרלו (התייחסות ל-50 ריצות)
mc = load('monte_carlo_results.mat');
fig1 = figure('Name', '1. Monte Carlo Analysis', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.05 0.5 0.4 0.4]);
tiledlayout(1,2);
nexttile;
imagesc(mc.r_grid/1000, mc.z_grid, mc.E_TL_mc);
title('Monte Carlo: E[TL]'); colorbar; colormap(gca, 'jet');
set(gca, 'YDir', 'reverse'); xlabel('Range (km)'); ylabel('Depth (m)');
nexttile;
imagesc(mc.r_grid/1000, mc.z_grid, mc.Var_TL_mc);
title('Monte Carlo: Var(TL)'); colorbar; colormap(gca, 'hot');
set(gca, 'YDir', 'reverse'); xlabel('Range (km)'); 
clim([0 prctile(mc.Var_TL_mc(:), 98, 'all')]); % חיתוך רעשי קצה לתצוגה נקייה

%% 2. תצוגת שיטת הנגזרת (Delta Method)
dm = load('delta_method_results.mat');
fig2 = figure('Name', '2. Delta Method Analysis', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.5 0.5 0.4 0.4]);
tiledlayout(1,2);
nexttile;
imagesc(dm.r_grid/1000, dm.z_grid, dm.E_TL_delta);
title('Delta Method: E[TL]'); colorbar; colormap(gca, 'jet');
set(gca, 'YDir', 'reverse'); xlabel('Range (km)'); ylabel('Depth (m)');
nexttile;
imagesc(dm.r_grid/1000, dm.z_grid, dm.Var_TL_delta);
title('Delta Method: Var(TL)'); colorbar; colormap(gca, 'hot');
set(gca, 'YDir', 'reverse'); xlabel('Range (km)');
clim([0 prctile(dm.Var_TL_delta(:), 98, 'all')]);

%% 3. השוואה: מונטה קרלו מול נגזרת
fig3 = figure('Name', '3. Comparison: MC vs Delta', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.05 0.05 0.4 0.4]);
tiledlayout(1,2);
% חישוב הפרשים (בהנחה שהרשתות זהות בזכות האינטרפולציה)
diff_E = mc.E_TL_mc - dm.E_TL_delta;
diff_V = mc.Var_TL_mc - dm.Var_TL_delta;
nexttile;
imagesc(mc.r_grid/1000, mc.z_grid, diff_E);
title('\Delta Expected Value (MC - Delta)'); colorbar; 
set(gca, 'YDir', 'reverse'); colormap(gca, parula);
limit_E = max(abs(diff_E(:))); clim([-limit_E limit_E]); % סקאלה סימטרית
nexttile;
imagesc(mc.r_grid/1000, mc.z_grid, diff_V);
title('\Delta Variance (MC - Delta)'); colorbar;
set(gca, 'YDir', 'reverse'); colormap(gca, parula);
limit_V = max(abs(diff_V(:))); clim([-limit_V limit_V]);

%% 4. הפעלת כלי אינטראקציה דינמי על כל הגרפים
allFigs = [fig0, fig1, fig2, fig3];
for f = 1:length(allFigs)
    dcm = datacursormode(allFigs(f));
    set(dcm, 'UpdateFcn', @customDatatipPresentation);
    datacursormode(allFigs(f), 'on');
end
msgbox('כל הגרפים מוכנים לפרזנטציה. השתמש בעכבר כדי לחקור נקודות בשדה.', 'Presentation Ready');

%% פונקציית דאטה-טיפ דינמית (מתאימה לכל גרף באופן אוטומטי)
function txt = customDatatipPresentation(~, event_obj)
    pos = get(event_obj, 'Position');
    img_obj = get(event_obj, 'Target');
    
    % חילוץ רשתות ה-X וה-Y ישירות מהגרף שעליו הקלקנו
    r_grid_km = get(img_obj, 'XData');
    z_grid = get(img_obj, 'YData');
    CData = get(img_obj, 'CData');
    
    % מציאת האינדקסים הקרובים ביותר ללחיצה
    [~, r_idx] = min(abs(r_grid_km - pos(1)));
    [~, z_idx] = min(abs(z_grid - pos(2)));
    current_val = CData(z_idx, r_idx);
    
    txt = {['Range: ', num2str(pos(1), '%.2f'), ' km'], ...
           ['Depth: ', num2str(pos(2), '%.1f'), ' m'], ...
           ['Value: ', num2str(current_val, '%.2f')]};
end