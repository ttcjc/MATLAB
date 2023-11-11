clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

figSave = false; % Save .fig File(s);

%%

mapDataA = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/SB_1.0L_120s_15Hz_01/T0067_T120000_F15_Norm.mat', 'mapData').mapData;
mapDataB = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/SB_1.0L_120s_15Hz_02/T0067_T120000_F15_Norm.mat', 'mapData').mapData;
mapDataC = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/SB_1.0L_120s_15Hz_03/T0067_T120000_F15_Norm.mat', 'mapData').mapData;


%%

[~, index] = max(mapDataB.density.mean);
y = mapDataA.positionGrid(index,2);
index = find(mapDataA.positionGrid(:,2) == mapDataA.positionGrid(index,2));


z = mapDataA.positionGrid(index,3);

rhoMeanA = mapDataA.density.mean(index,:);
rhoMeanB = mapDataB.density.mean(index,:);
rhoMeanC = mapDataC.density.mean(index,:);

rhoRMSa = mapDataA.density.RMS(index,:);
rhoRMSb = mapDataB.density.RMS(index,:);
rhoRMSc = mapDataC.density.RMS(index,:);


%%

% Initialise Figure
fig = fig + 1;
figName = ['Mean_Spray_Density_Profile_y_', num2str(y)];
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 20, 'layer', 'top');

% Plot Mean Profiles
plot(z, rhoMeanA, 'color', graphColours(1), 'lineWidth', 2);
plot(z, rhoMeanB, 'color', graphColours(2), 'lineWidth', 2);
plot(z, rhoMeanC, 'color', graphColours(3), 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([0; 0.5]);
ylim([0; 1.05]);
tickData = (0.1:0.1:0.4);
xticks(tickData);
tickData = (0.21:0.21:0.84);
yticks(tickData);
xtickformat('%.1f');
ytickformat('%.2f');
xlabel({'{$z_{\ell}$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$\bar{\varrho_{n}}$}'}, 'interpreter', 'latex');
tightInset = get(gca, 'TightInset');
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
figName = ['RMS_Spray_Density_Profile_y_', num2str(z(1))];
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.775], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 20, 'layer', 'top');

% Plot RMS Profiles
plot(z, rhoRMSa, 'color', graphColours(1), 'lineWidth', 2);
plot(z, rhoRMSb, 'color', graphColours(2), 'lineWidth', 2);
plot(z, rhoRMSc, 'color', graphColours(3), 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([0; 0.5]);
ylim([0; 0.7]);
tickData = (0.1:0.1:0.4);
xticks(tickData);
tickData = (0.14:0.14:0.56);
yticks(tickData);
xtickformat('%.1f');
ytickformat('%.2f');
xlabel({'{$z_{\ell}$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$RMS(\varrho_{n})$}'}, 'interpreter', 'latex');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');
