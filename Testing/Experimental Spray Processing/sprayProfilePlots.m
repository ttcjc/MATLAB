run preamble;

normLength = 1.044;
refValue = 0.0052166;

figSave = false; % Save .fig File(s);

%%

caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.0L_120s_15Hz_01/T0067_T120000_F15.mat';
caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.0L_120s_15Hz_02/T0067_T120000_F15.mat';
caseC = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.0L_120s_15Hz_03/T0067_T120000_F15.mat';
caseD = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.0L_120s_15Hz_04/T0067_T120000_F15.mat';
% caseE = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.0L_600s_03Hz_01/T0333_T600000_F3.mat';    

% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_1.0L_120s_15Hz_01/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_1.0L_120s_15Hz_02/T0067_T120000_F15.mat';
% caseC = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_1.0L_120s_15Hz_03/T0067_T120000_F15.mat';
% caseD = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_1.0L_120s_15Hz_04/T0067_T120000_F15.mat';
% caseE = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_1.0L_600s_03Hz_01/T0333_T600000_F3.mat';    

% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_1.0L_120s_15Hz_01/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_1.0L_120s_15Hz_02/T0067_T120000_F15.mat';
% caseC = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_1.0L_120s_15Hz_03/T0067_T120000_F15.mat';
% caseD = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_1.0L_120s_15Hz_04/T0067_T120000_F15.mat';
% caseE = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_1.0L_600s_03Hz_01/T0333_T600000_F3.mat';    

% Load Data
mapDataA = load(caseA, 'mapData').mapData;
mapDataB = load(caseB, 'mapData').mapData;
mapDataC = load(caseC, 'mapData').mapData;
mapDataD = load(caseD, 'mapData').mapData;
% mapDataE = load(caseE, 'mapData').mapData;

% sprayMaps = fieldnames(R2R);


%%

index = 4672; % [~, index] = max(mapDataB.density.mean);
y = mapDataA.positionGrid(index,2) / normLength;
index = find(mapDataA.positionGrid(:,2) == mapDataA.positionGrid(index,2));
z = mapDataA.positionGrid(index,3) / normLength;

rhoMeanA = mapDataA.density.mean(index,:) / refValue;
rhoMeanB = mapDataB.density.mean(index,:) / refValue;
rhoMeanC = mapDataC.density.mean(index,:) / refValue;
rhoMeanD = mapDataD.density.mean(index,:) / refValue;

rhoRMSa = mapDataA.density.RMS(index,:) / refValue;
rhoRMSb = mapDataB.density.RMS(index,:) / refValue;
rhoRMSc = mapDataC.density.RMS(index,:) / refValue;
rhoRMSd = mapDataD.density.RMS(index,:) / refValue;


%%

% Initialise Figure
fig = fig + 1;
figName = ['Run_2_Run_Variation_Mean_Y=', num2str(y)];
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot Mean Profiles
plot(z, rhoMeanA, 'color', graphColours(1), 'lineWidth', 2);
plot(z, rhoMeanB, 'color', graphColours(2), 'lineWidth', 2);
plot(z, rhoMeanC, 'color', graphColours(3), 'lineWidth', 2);
plot(z, rhoMeanD, 'color', graphColours(4), 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([0; 0.5]);
ylim([0; 1.1]);
tickData = (0.1:0.1:0.4);
xticks(tickData);
tickData = (0.22:0.22:0.88);
yticks(tickData);
xtickformat('%.1f');
ytickformat('%.2f');
xlabel({'{$z_{\ell}$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$\overline{\varrho_{_{n}}}$}'}, 'interpreter', 'latex');
tightInset = get(gca, 'TightInset');
legend({'Run 1', ...
        'Run 2', ...
        'Run 3', ...
        'Run 4'}, ...
       'location', 'northEast', 'orientation', 'vertical', 'interpreter', 'latex', ...
       'fontSize', 18, 'box', 'off');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');


%%

% Initialise Figure
fig = fig + 1;
figName = ['Run_2_Run_Variation_RMS_Y=', num2str(y)];
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot RMS Profiles
plot(z, rhoRMSa, 'color', graphColours(1), 'lineWidth', 2);
plot(z, rhoRMSb, 'color', graphColours(2), 'lineWidth', 2);
plot(z, rhoRMSc, 'color', graphColours(3), 'lineWidth', 2);
plot(z, rhoRMSd, 'color', graphColours(4), 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([0; 0.5]);
ylim([0; 0.75]);
tickData = (0.1:0.1:0.4);
xticks(tickData);
tickData = (0.15:0.15:0.6);
yticks(tickData);
xtickformat('%.1f');
ytickformat('%.2f');
xlabel({'{$z_{\ell}$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$\mathrm{RMS}(\varrho_{_{n}})$}'}, 'interpreter', 'latex');
tightInset = get(gca, 'TightInset');
legend({'Run 1', ...
        'Run 2', ...
        'Run 3', ...
        'Run 4'}, ...
       'location', 'northEast', 'orientation', 'vertical', 'interpreter', 'latex', ...
       'fontSize', 18, 'box', 'off');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');
