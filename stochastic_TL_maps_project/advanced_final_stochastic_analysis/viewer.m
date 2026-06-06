%% viewer.m  –  Interactive TL UQ Dashboard  (2026-05-26)
%
% FEATURES
% ────────
%   • Scenario selector  (Baseline / Deep Water / Shallow Water / Upslope / Downslope)
%   • Subset selector    (7 parameter subsets)
%   • Display mode       (Mean TL | Variance | Side-by-Side)
%   • Probability map    (continuous P(shadow) with threshold contours)
%   • Method selector    (Delta Cheb / MC Cheb / MC Empirical / MC LN3)
%   • Threshold overlays (70% / 80% / 90% / 95%)
%   • DataTip            (hover shows Mean + Variance + Std at any point)

try, cd(fileparts(mfilename('fullpath'))); catch; end
addpath(genpath('Functions'));

%% ── Data registry ───────────────────────────────────────────────────────
SCENARIOS = { ...
    'Baseline (const_35, 50km)',   'Delta_Method/results', 'Monte_Carlo/results'; ...
    'Deep Water (2500m, 50km)',    'Scenarios/deep_water/results',  'Scenarios/deep_water/results'; ...
    'Shallow Water (50m, 5km)',    'Scenarios/shallow_water/results','Scenarios/shallow_water/results'; ...
    'Upslope (550→50m, 10km)',     'Scenarios/upslope/results',     'Scenarios/upslope/results'; ...
    'Downslope (50→550m, 10km)',   'Scenarios/downslope/results',   'Scenarios/downslope/results'; ...
};

subset_names  = {'zS','freq','svp','zS_freq','zS_svp','freq_svp','zS_freq_svp'};
subset_labels = {'z_S only','Freq only','SVP only','z_S+Freq','z_S+SVP','Freq+SVP','z_S+Freq+SVP'};

PROB_METHODS = {'Delta Chebyshev','MC Chebyshev','MC Empirical','MC LN3'};
THRESHOLDS   = [0.70 0.80 0.90 0.95];
FOM          = 100;

%% ── Build UI ────────────────────────────────────────────────────────────
fig = uifigure('Name','Stochastic Bellhop Viewer', ...
               'Position',[50 50 1500 800]);

% Left control panel
cpanel = uipanel(fig,'Position',[5 5 190 790],'Title','Controls');
y0=730;
function h=lbl(txt,y), h=uilabel(cpanel,'Text',txt,'Position',[8 y 174 20],'FontWeight','bold'); end
function h=dd(items,val,y), h=uidropdown(cpanel,'Items',items,'Value',val,'Position',[8 y 174 28]); end

lbl('Scenario',y0);           dd_scn  = dd(SCENARIOS(:,1)', SCENARIOS{1,1},  y0-32);
lbl('Subset',y0-70);          dd_sub  = dd(subset_labels,   subset_labels{7},y0-102);
lbl('Display Mode',y0-140);
dd_mode = uidropdown(cpanel,'Items',{'Mean TL','Variance','Side-by-Side'}, ...
    'Value','Mean TL','Position',[8 y0-172 174 28]);
lbl('Prob Method',y0-210);    dd_prob = dd(PROB_METHODS,     'MC Empirical',  y0-242);
lbl('Show Prob Map',y0-280);
chk_prob = uicheckbox(cpanel,'Text','enable','Value',false,'Position',[8 y0-302 174 22]);
uibutton(cpanel,'Text','Update','Position',[8 y0-345 174 34], ...
    'ButtonPushedFcn',@(~,~)updatePlot(),'BackgroundColor',[0.3 0.6 1],'FontWeight','bold');

% Right axes area
AX_W=600; AX_H=500; AX_Y=240;
ax1=uiaxes(fig,'Position',[205  AX_Y  AX_W  AX_H]);
ax2=uiaxes(fig,'Position',[830  AX_Y  AX_W  AX_H],'Visible','off');
ax3=uiaxes(fig,'Position',[205  20    1280  200],  'Visible','off');

%% ── State ───────────────────────────────────────────────────────────────
D=[]; M=[]; cache_scn=-1; cache_sub=-1;

%% ── Load data ────────────────────────────────────────────────────────────
    function loadData()
        si=find(strcmp(SCENARIOS(:,1),dd_scn.Value)); if isempty(si),si=1;end
        bi=find(strcmp(subset_labels,dd_sub.Value));   if isempty(bi),bi=7;end
        if si==cache_scn&&bi==cache_sub, return; end
        sn=subset_names{bi};
        df=fullfile(SCENARIOS{si,2},sprintf('delta_%s.mat',sn));
        mf=fullfile(SCENARIOS{si,3},sprintf('MC_%s.mat',sn));
        D=[]; M=[];
        if isfile(df), D=load(df); end
        if isfile(mf), M=load(mf); end
        cache_scn=si; cache_sub=bi;
        if isempty(D)&&isempty(M)
            uialert(fig,sprintf('No results found.\nScenario: %s\nSubset: %s\n\nRun the analysis scripts first.', ...
                SCENARIOS{si,1},subset_labels{bi}),'Data Missing');
        end
    end

%% ── Update plot ──────────────────────────────────────────────────────────
    function updatePlot()
        loadData();
        if isempty(D)&&isempty(M), return; end

        % Pick canonical data source (prefer MC for mean/var, Delta for TL_expected)
        if ~isempty(D)
            r_km=D.r_km; z_m=D.z_m;
            mean_map=D.TL_expected; var_map=D.Var_TL;
            % Override with MC if available (more accurate mean/var)
            if ~isempty(M)&&isfield(M,'MC_EX')
                mean_map=M.MC_EX; var_map=M.MC_Var;
            end
        else
            r_km=M.r_km; z_m=M.z_m;
            mean_map=M.MC_EX; var_map=M.MC_Var;
        end
        TL_contour=movmean(movmean(mean_map,20,2),20,1);
        sub_lbl=dd_sub.Value;

        switch dd_mode.Value
            case 'Mean TL'
                ax2.Visible='off'; cla(ax2);
                plotMean(ax1,r_km,z_m,mean_map,var_map,TL_contour,sub_lbl);
            case 'Variance'
                ax2.Visible='off'; cla(ax2);
                plotVar(ax1,r_km,z_m,var_map,sub_lbl);
            case 'Side-by-Side'
                ax2.Visible='on';
                plotMean(ax1,r_km,z_m,mean_map,var_map,TL_contour,['Mean | ' sub_lbl]);
                plotVar( ax2,r_km,z_m,var_map,['Var | ' sub_lbl]);
        end

        % Probability pane
        if chk_prob.Value
            ax3.Visible='on';
            pm=getProbMap(); r2=r_km; z2=z_m;
            if ~isempty(pm)
                cla(ax3);
                shadowThresholdMaps(ax3,r2,z2,pm, ...
                    sprintf('P(TL>%ddB) | %s | %s',FOM,dd_prob.Value,sub_lbl), ...
                    THRESHOLDS);
            end
        else
            ax3.Visible='off'; cla(ax3);
        end
    end

%% ── Plot functions ───────────────────────────────────────────────────────
    function plotMean(ax,r_km,z_m,mean_map,var_map,TL_c,ttl)
        cla(ax);
        pcolor(ax,r_km,z_m,mean_map); shading(ax,'interp');
        set(ax,'YDir','reverse'); colormap(ax,jet); clim(ax,[50 150]); colorbar(ax);
        hold(ax,'on');
        contour(ax,r_km,z_m,TL_c,[FOM FOM],'w-','LineWidth',0.8);
        hold(ax,'off');
        xlabel(ax,'Range (km)'); ylabel(ax,'Depth (m)'); title(ax,ttl,'FontSize',9);
        % DataTip with mean + variance
        dcm=datacursormode(fig);
        set(dcm,'UpdateFcn',@(~,e)meanVarTip(e,r_km,z_m,mean_map,var_map));
    end

    function plotVar(ax,r_km,z_m,var_map,ttl)
        cla(ax);
        pcolor(ax,r_km,z_m,var_map); shading(ax,'interp');
        set(ax,'YDir','reverse'); colormap(ax,hot); colorbar(ax);
        xlabel(ax,'Range (km)'); ylabel(ax,'Depth (m)'); title(ax,ttl,'FontSize',9);
    end

    function pm=getProbMap()
        pm=[];
        switch dd_prob.Value
            case 'Delta Chebyshev', if ~isempty(D)&&isfield(D,'Cheb_lb_s'),pm=D.Cheb_lb_s;end
            case 'MC Chebyshev',    if ~isempty(M)&&isfield(M,'Cheb_lb'),  pm=M.Cheb_lb;  end
            case 'MC Empirical',    if ~isempty(M)&&isfield(M,'MC_PrFOM'), pm=M.MC_PrFOM; end
            case 'MC LN3',          if ~isempty(M)&&isfield(M,'LN3_prob'), pm=M.LN3_prob; end
        end
    end

    function txt=meanVarTip(evt,r_km,z_m,mean_map,var_map)
        pos=get(evt,'Position');
        [~,ri]=min(abs(r_km-pos(1))); [~,zi]=min(abs(z_m-pos(2)));
        ri=max(1,min(ri,numel(r_km))); zi=max(1,min(zi,numel(z_m)));
        txt={sprintf('Range: %.2f km',  r_km(ri)), ...
             sprintf('Depth: %.1f m',   z_m(zi)), ...
             sprintf('Mean TL: %.3f dB',mean_map(zi,ri)), ...
             sprintf('Var[TL]: %.4f dB²',var_map(zi,ri)), ...
             sprintf('Std[TL]: %.4f dB',sqrt(max(0,var_map(zi,ri))))};
    end

% Initial render
updatePlot();
