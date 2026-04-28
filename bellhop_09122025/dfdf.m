clear; clc; close all;

lambda = 1.55e-6; 
delta_H = pi;     
delta_L = pi;
N_max = 500;      
N_range = 1:N_max;

ratios = [1.1, 1.01];
colors = {'b', 'm'}; 
labels = {'n_h/n_l = 1.1', 'n_h/n_l = 1.01'};
guess_points = [20, 200]; 

figure('Color', 'w', 'Position', [100, 100, 900, 600]);
hold on; grid on;

for k = 1:length(ratios)
    ratio = ratios(k); 
    
    r = (ratio - 1) / (ratio + 1); 
    
    A_val = exp(-1i*delta_L/2) * (exp(-1i*delta_H/2) - r^2 * exp(1i*delta_H/2)); 
    B_val = exp(-1i*delta_L/2) * (-r * exp(-1i*delta_H/2) + r * exp(1i*delta_H/2)); %
    C_val = exp(1i*delta_L/2)  * (-r * exp(1i*delta_H/2) + r * exp(-1i*delta_H/2)); % 
    D_val = exp(1i*delta_L/2)  * (exp(1i*delta_H/2) - r^2 * exp(-1i*delta_H/2)); % 
    
    M_base = (1 / (1 - r^2)) * [A_val, B_val; C_val, D_val]; 
    
    R = zeros(size(N_range));
    for i = 1:length(N_range)
        n = N_range(i);
        Mn = M_base^n;
        C_n = Mn(2,1);
        D_n = Mn(2,2); 
        R(i) = abs(C_n / D_n)^2;
    end
    
    plot(N_range, R, 'Color', colors{k}, 'LineWidth', 2, 'DisplayName', labels{k});
    
    x_guess = guess_points(k);
    y_guess = R(N_range == x_guess);
    plot(x_guess, y_guess, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
    text(x_guess, y_guess + 0.03, sprintf('guess for ratio %.2f (N=%d)', ratio, x_guess), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 9);
    
    [~, idx] = min(abs(R - 0.99));
    n_99 = N_range(idx);
    r_99 = R(idx);
    
    plot(n_99, r_99, 'kx', 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
    text(n_99, r_99 - 0.05, sprintf('R=0.99 at N=%d', n_99), ...
        'Color', 'k', 'FontSize', 9, 'VerticalAlignment', 'top');
end

title('reflection vs N', 'FontSize', 14);
xlabel('N (Number of periods)', 'FontSize', 12);
ylabel('(C/D)^2 (Reflection)', 'FontSize', 12);
legend('Location', 'southeast');
ylim([0 1.1]);
hold off;