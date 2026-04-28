% ========================================================================
% Tight-Binding Dispersion Relation for 1D Benzene Chain
% Plotting the First Brillouin Zone
% ========================================================================
clear; clc; close all;

%% 1. Physical Parameters
% Using typical values for graphene-like systems for realism
t = 2.7;          % Hopping integral [eV]
alpha = 0;        % On-site energy [eV] (Reference level, sets Ef approx at 0)
a = 1.42;         % Bond length between Carbons [Angstrom]

% Lattice constant in real space
R = a * sqrt(3);  % [Angstrom]

%% 2. Defining the First Brillouin Zone (BZ)
% The BZ boundaries are at +/- pi/R
k_boundary = pi / R;

% Create k-vector spanning exactly the first BZ
num_points = 1000;
k = linspace(-k_boundary, k_boundary, num_points);

%% 3. Calculating Energy Bands
% Recall the formula derived: E = -alpha +/- (t/4)*( +/-1 +/- S )
% Where S = sqrt(1 + 16*cos^2(kR/2))

% Helper term S(k)
S_k = sqrt(1 + 16 * (cos(k * R / 2)).^2);

% The four energy bands (ordered from lowest to highest)
% Band 1 (Lowest)
E1 = -alpha - (t/4) * (S_k + 1);
% Band 2 (Valence Band - Highest occupied at 0K)
E2 = -alpha + (t/4) * (1 - S_k);
% Band 3 (Conduction Band - Lowest unoccupied at 0K)
E3 = -alpha + (t/4) * (S_k - 1);
% Band 4 (Highest)
E4 = -alpha + (t/4) * (S_k + 1);

%% 4. Plotting
figure('Name', '1D Benzene Chain Dispersion', 'Color', 'w', 'Position', [100, 100, 800, 600]);
hold on; box on;

% Plot the bands with distinct styles
plot(k, E1, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Band 1 (Lowest)');
plot(k, E2, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Band 2 (Valence)');
plot(k, E3, 'g-', 'LineWidth', 2.5, 'DisplayName', 'Band 3 (Conduction)');
plot(k, E4, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Band 4 (Highest)');

% Mark the Fermi Energy level (approximate for this symmetric case)
yline(-alpha, 'k:', 'Fermi Energy (E_F)', 'LabelVerticalAlignment', 'middle', 'LineWidth', 1);

% --- Marking the Brillouin Zone Boundaries ---
xline(-k_boundary, 'k-.', {'-pi/R', 'BZ Boundary'}, 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.2);
xline(k_boundary, 'k-.', {'+pi/R', 'BZ Boundary'}, 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.2);

% Graph decorations
title({'Electronic Dispersion Relation', '1D Benzene Chain (First Brillouin Zone)'}, 'FontSize', 14);
xlabel('Wavevector k [1/\AA]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Energy E(k) [eV]', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'bestoutside', 'FontSize', 11);
grid on;

% Set x-axis limits exactly to the BZ
xlim([-k_boundary, k_boundary]);

% Add a slight padding to y-axis for better view
ylim([min(E1)*1.1, max(E4)*1.1]);

set(gca, 'FontSize', 10, 'LineWidth', 1.2);
hold off;