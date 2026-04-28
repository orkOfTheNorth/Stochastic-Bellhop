%% Phase Diagram: Metal-Insulator Transition
clear; clc; clf;

% 1. הגדרת הפרמטר הקבוע
% זהו הפרש האנרגיות באתר בין אורביטל B ל-A
delta_eps = 4; % ערך לדוגמה (epsilon^B - epsilon^A)

% 2. יצירת רשת של ערכי t_AA ו-t_BB
grid_points = 500;
t_range = linspace(0, 1.5 * delta_eps, grid_points);
[T_AA, T_BB] = meshgrid(t_range, t_range);

% 3. חישוב מצב הפאזה עבור כל נקודה ברשת
% התנאי למתכת: T_AA + T_BB > delta_eps
% המטריצה 'PhaseState' תכיל 0 עבור מבודד ו-1 עבור מתכת
PhaseState = (T_AA + T_BB) > delta_eps;

% 4. יצירת השרטוט
figure('Color', 'w', 'Name', 'Phase Diagram');

% שימוש ב-imagesc כדי לצבוע את האזורים
% אנו משתמשים ב-PhaseState כדי לקבוע את הצבע
imagesc(t_range, t_range, PhaseState);
hold on;

% הגדרת כיוון ציר Y (שלא יהיה הפוך כמו בתמונה רגילה)
set(gca, 'YDir', 'normal');

% הגדרת מפת צבעים: כחול למבודד (0), אדום למתכת (1)
colormap([0.7 0.8 1; 1 0.7 0.7]); 

% 5. הוספת קו הגבול
% קו הגבול הוא הישר: t_BB = delta_eps - t_AA
plot([0, delta_eps], [delta_eps, 0], 'k-', 'LineWidth', 3);

% 6. הוספת כותרות ותוויות
title(['דיאגרמת פאזה: מתכת-מבודד (\Delta\epsilon = ', num2str(delta_eps), ')'], ...
    'FontSize', 14);
xlabel('מקדם דילוג t_{AA}', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('מקדם דילוג t_{BB}', 'FontSize', 12, 'FontWeight', 'bold');

% הוספת טקסט לציון האזורים
text(delta_eps/3, delta_eps/3, 'INSULATOR (מבודד)', ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
text(delta_eps*1.1, delta_eps*1.1, 'METAL (מתכת)', ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');

axis square; % הופך את הצירים לריבועיים
grid on;
set(gca, 'Layer', 'top'); % שם את הרשת וקו הגבול מעל הצבע