clear; clc; close all;

lambda0 = 1.55; 
nL = 1.45;
lambda_vec = linspace(1.45, 1.65, 2000);
ratios = [1.01, 1.1]; 
N_vals = [301, 32]; 

figure('Name', 'DBR Analysis', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

for i = 1:2
    nH = nL * ratios(i);
    N = N_vals(i);
    r = (nL - nH) / (nL + nH);
    
    dL = lambda0 / (4 * nL);
    dH = lambda0 / (4 * nH);
    
    R_results = zeros(size(lambda_vec));
    Phase_results = zeros(size(lambda_vec));
    
    for k = 1:length(lambda_vec)
        lam = lambda_vec(k);
        deltaL = (4 * pi * nL * dL) / lam;
        deltaH = (4 * pi * nH * dH) / lam;
        
        A = exp(-1j*deltaL/2) * (exp(-1j*deltaH/2) - r^2 * exp(1j*deltaH/2));
        B = exp(-1j*deltaL/2) * (-r*exp(-1j*deltaH/2) + r*exp(1j*deltaH/2));
        C = exp(1j*deltaL/2) * (-r*exp(1j*deltaH/2) + r*exp(-1j*deltaH/2));
        D = exp(1j*deltaL/2) * (exp(1j*deltaH/2) - r^2 * exp(-1j*deltaH/2));
        
        M_single = (1 / (1 - r^2)) * [A, B; C, D];
        M_total = M_single^N;
        
        M_exit = (1/(1+r)) * [1, r; r, 1];
        M_sys = M_total * M_exit;
        
        rN = -M_sys(2,1) / M_sys(2,2);
        
        R_results(k) = abs(rN)^2;
        Phase_results(k) = angle(rN);
    end
    
    subplot(2, 2, 2*i - 1);
    plot(lambda_vec, R_results, 'LineWidth', 2); hold on;
    yline(0.99, 'r--', '99%', 'LabelHorizontalAlignment', 'left');
    grid on;
    title(['Reflectivity (R): n_H/n_L = ', num2str(ratios(i)), ', N = ', num2str(N)]);
    xlabel('\lambda [\mum]'); ylabel('Reflectivity');
    ylim([0 1.05]); 
    xlim([1.45 1.65]);

    subplot(2, 2, 2*i);
    plot(lambda_vec, unwrap(Phase_results), 'Color', [0.85 0.325 0.098], 'LineWidth', 2);
    grid on;
    title(['Phase (Unwrapped): n_H/n_L = ', num2str(ratios(i))]);
    xlabel('\lambda [\mum]'); ylabel('Phase [rad]');
    xlim([1.45 1.65]);
end

sgtitle('DBR Performance Summary - Peak Reflectivity Analysis');