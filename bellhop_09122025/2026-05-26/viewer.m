%% viewer.m  –  Master UQ Results Dashboard  (2026-05-26)
%
% Interactive GUI to browse all saved uncertainty maps.
% Select method, subset, and map type from the control panel.
%
% Workflow:
%   1. Run  Delta_Method/run_delta.m
%   2. Run  Monte_Carlo/run_MC.m
%   3. Run  Comparison/run_comparison.m
%   4. Open this viewer to explore all outputs interactively.

clear; close all; clc;
try, cd(fileparts(mfilename('fullpath'))); catch; end

FOM = 100;

subset_names  = {'zS','freq','temp','zS_freq','zS_temp','freq_temp','zS_freq_temp'};
subset_labels = {'z_S only','Freq only','Temp only', ...
                 'z_S + Freq','z_S + Temp','Freq + Temp','z_S + Freq + Temp'};

%% ── UI FIGURE ────────────────────────────────────────────────────────────────
ui = uifigure('Name','UQ Results Viewer – 2026-05-26','Position',[50 50 1100 780]);

% Title
uilabel(ui,'Position',[20 740 1060 30],'Text', ...
    'UQ Results Viewer  |  freq±1%  zS±1%  Temp±1°C  |  FOM = 100 dB', ...
    'FontSize',14,'FontWeight','bold','HorizontalAlignment','center');

% Method selector
uilabel(ui,'Position',[20 700 80 22],'Text','Method:');
dd_method = uidropdown(ui,'Position',[100 700 150 25], ...
    'Items',{'Delta Method','Monte Carlo','Comparison'}, ...
    'Value','Delta Method');

% Subset selector
uilabel(ui,'Position',[270 700 55 22],'Text','Subset:');
dd_subset = uidropdown(ui,'Position',[330 700 200 25], ...
    'Items', subset_labels, 'Value', subset_labels{1});

% Map type selector
uilabel(ui,'Position',[550 700 65 22],'Text','Map type:');
dd_maptype = uidropdown(ui,'Position',[620 700 200 25], ...
    'Items',{'E[TL] (Expected)','Var[TL] (Variance)','Chebyshev Shadow %', ...
             'MC Empirical P(TL>FOM)','Shadow Zone (binary 95%)'}, ...
    'Value','E[TL] (Expected)');

% Load & Display button
uibutton(ui,'push','Position',[840 695 120 32],'Text','Load & Display', ...
    'ButtonPushedFcn', @(~,~) loadAndDisplay(dd_method.Value, dd_subset.Value, ...
        dd_maptype.Value, subset_names, subset_labels, FOM));

% Open MC interactive button
uibutton(ui,'push','Position',[975 695 105 32],'Text','MC Interactive', ...
    'FontColor',[0.1 0.4 0.8], ...
    'ButtonPushedFcn', @(~,~) runMCViewer());

% Axes for main map
ax_main = uiaxes(ui,'Position',[20 60 760 620]);

% Info panel
uilabel(ui,'Position',[800 640 280 22],'Text','Summary Info:','FontWeight','bold');
txt_info = uitextarea(ui,'Position',[800 60 280 570],'Editable','off', ...
    'Value',{'Select a result and click [Load & Display]', '', ...
             'Method order:', '  1. Run Delta_Method/run_delta.m', ...
             '  2. Run Monte_Carlo/run_MC.m', ...
             '  3. Run Comparison/run_comparison.m', '', ...
             'Then browse results here.'});

%% ── CALLBACK: Load and display selected result ───────────────────────────────
    function loadAndDisplay(method, subset_lbl, maptype, snames, slabels, fom)
        % Find subset index
        s_idx = find(strcmp(slabels, subset_lbl));
        if isempty(s_idx), return; end
        sn = snames{s_idx};

        try
            switch method
                case 'Delta Method'
                    data = load(fullfile('Delta_Method','results', ...
                                         sprintf('delta_%s.mat',sn)));
                    switch maptype
                        case 'E[TL] (Expected)'
                            plotMap(data.TL_expected, data.r_km, data.z_m, jet, ...
                                [50 150], sprintf('Delta E[TL] | %s',subset_lbl), fom);
                        case 'Var[TL] (Variance)'
                            plotMap(data.Var_TL, data.r_km, data.z_m, hot, ...
                                [], sprintf('Delta Var[TL] | %s',subset_lbl), fom);
                        case 'Chebyshev Shadow %'
                            plotMap(data.Cheb_lb_s*100, data.r_km, data.z_m, parula, ...
                                [0 100], sprintf('Delta Chebyshev P(TL>%ddB) %% | %s',fom,subset_lbl), fom);
                        case 'Shadow Zone (binary 95%)'
                            plotMap(double(data.Cheb_shadow95), data.r_km, data.z_m, ...
                                [0.85 0.85 0.85; 0.2 0.7 0.2], [0 1], ...
                                sprintf('Delta Shadow Zone 95%% | %s',subset_lbl), fom);
                        otherwise
                            txt_info.Value = {'Map type not available for Delta method.'};
                            return;
                    end
                    txt_info.Value = formatInfo('Delta', subset_lbl, sn, data, fom, maptype);

                case 'Monte Carlo'
                    data = load(fullfile('Monte_Carlo','results', ...
                                         sprintf('MC_%s.mat',sn)));
                    switch maptype
                        case 'E[TL] (Expected)'
                            plotMap(data.MC_EX, data.r_km, data.z_m, jet, ...
                                [50 150], sprintf('MC E[TL] | N=%d | %s',data.N,subset_lbl), fom);
                        case 'Var[TL] (Variance)'
                            plotMap(data.MC_Var, data.r_km, data.z_m, hot, ...
                                [], sprintf('MC Var[TL] | N=%d | %s',data.N,subset_lbl), fom);
                        case 'Chebyshev Shadow %'
                            plotMap(data.Cheb_lb*100, data.r_km, data.z_m, parula, ...
                                [0 100], sprintf('MC Chebyshev P(TL>%ddB) %% | %s',fom,subset_lbl), fom);
                        case 'MC Empirical P(TL>FOM)'
                            plotMap(data.MC_PrFOM*100, data.r_km, data.z_m, parula, ...
                                [0 100], sprintf('MC Empirical P(TL>%ddB) %% | N=%d | %s',fom,data.N,subset_lbl), fom);
                        case 'Shadow Zone (binary 95%)'
                            plotMap(double(data.MC_PrFOM>=0.95), data.r_km, data.z_m, ...
                                [0.85 0.85 0.85; 0.2 0.7 0.2], [0 1], ...
                                sprintf('MC Shadow Zone 95%% (empirical) | %s',subset_lbl), fom);
                    end
                    txt_info.Value = formatInfo('MC', subset_lbl, sn, data, fom, maptype);

                case 'Comparison'
                    data = load(fullfile('Comparison','results', ...
                                         sprintf('compare_%s.mat',sn)));
                    switch maptype
                        case 'E[TL] (Expected)'
                            plotMap(data.diff_EX, data.r_km, data.z_m, ...
                                redblue_local, [], ...
                                sprintf('MC−Delta E[TL] Diff | %s',subset_lbl), fom);
                        case 'Var[TL] (Variance)'
                            plotMap(data.diff_Var, data.r_km, data.z_m, ...
                                redblue_local, [], ...
                                sprintf('MC−Delta Var[TL] Diff | %s',subset_lbl), fom);
                        case 'Shadow Zone (binary 95%)'
                            shd_map = double(data.shd_delta_cheb) + 2*double(data.shd_mc_emp);
                            plotMap(shd_map, data.r_km, data.z_m, ...
                                [0.9 0.9 0.9; 0.2 0.5 0.9; 0.9 0.3 0.2; 0.2 0.7 0.2], [0 3], ...
                                sprintf('Shadow Zones | %s | Grey=none Blue=DeltaOnly Red=MConly Green=both',subset_lbl), fom);
                        otherwise
                            txt_info.Value = {'Use E[TL] or Var[TL] or Shadow for comparison.'};
                            return;
                    end
                    txt_info.Value = formatInfo('Compare', subset_lbl, sn, data, fom, maptype);
            end
        catch ME
            txt_info.Value = {['Error: ' ME.message], '', ...
                'Make sure you have run the method scripts first.'};
        end
    end

    function plotMap(data_map, r_km, z_m, cmap, clims, ttl, fom)
        fig = figure('Name',ttl,'Position',[200 100 900 620]);
        ax  = axes('Parent',fig);
        pcolor(ax, r_km, z_m, data_map); shading(ax,'interp');
        set(ax,'YDir','reverse');
        colormap(ax, cmap); colorbar(ax);
        if ~isempty(clims), clim(ax, clims); end
        xlabel(ax,'Range (km)'); ylabel(ax,'Depth (m)');
        title(ax, sprintf('%s\nFOM=%ddB  f=10kHz  zS=5m  freq±1%%  zS±1%%  Temp±1°C', ttl, fom));
    end

    function lines = formatInfo(method, subset_lbl, sn, data, fom, maptype)
        lines = {sprintf('Method:  %s', method), ...
                 sprintf('Subset:  %s', subset_lbl), ...
                 sprintf('Map:     %s', maptype), ...
                 sprintf('FOM:     %d dB', fom), ...
                 '', 'Statistics (spatial):'};
        if isfield(data,'Var_TL')
            lines{end+1} = sprintf('Max Var[TL]: %.4f', max(data.Var_TL(:)));
            lines{end+1} = sprintf('Mean Var[TL]: %.4f', mean(data.Var_TL(:)));
        elseif isfield(data,'MC_Var')
            lines{end+1} = sprintf('Max MC Var:  %.4f', max(data.MC_Var(:)));
            lines{end+1} = sprintf('Mean MC Var: %.4f', mean(data.MC_Var(:)));
        end
        if isfield(data,'Cheb_lb_s')
            lb = data.Cheb_lb_s;
            lines{end+1} = sprintf('Shadow area (Cheb 95%%): %.1f%%', ...
                100*mean(lb(:)>=0.95));
        elseif isfield(data,'Cheb_lb') && isnumeric(data.Cheb_lb)
            lb = data.Cheb_lb;
            lines{end+1} = sprintf('Shadow area (Cheb 95%%): %.1f%%', ...
                100*mean(lb(:)>=0.95));
        end
        if isfield(data,'MC_PrFOM')
            lines{end+1} = sprintf('Shadow area (MC emp 95%%): %.1f%%', ...
                100*mean(data.MC_PrFOM(:)>=0.95));
        end
        if isfield(data,'N')
            lines{end+1} = sprintf('MC runs N: %d', data.N);
        end
    end

    function runMCViewer()
        run('Monte_Carlo/show_MC.m');
    end

    function cmap = redblue_local()
        n = 256;
        r1 = [linspace(0.2,1,n/2); linspace(0.2,1,n/2); ones(1,n/2)]';
        r2 = [ones(1,n/2); linspace(1,0.2,n/2); linspace(1,0.2,n/2)]';
        cmap = [r1; r2];
    end
