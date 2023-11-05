clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

%%%

mapDataA = load('/mnt/Processing/Data/Experimental/MATLAB/planarContaminantMap/SB_1.0L_120s_15Hz_03/T0067_T120000_F15_Norm.mat', 'mapData').mapData;
mapDataB = load('/mnt/Processing/Data/Experimental/MATLAB/planarContaminantMap/SB_1.0L_120s_15Hz_02/T0067_T120000_F15_Norm.mat', 'mapData').mapData;
mapDataC = load('/mnt/Processing/Data/Experimental/MATLAB/planarContaminantMap/SB_1.0L_120s_15Hz_01/T0067_T120000_F15_Norm.mat', 'mapData').mapData;

%%

desiredPos = (-0.167 / 1.044);
[~, index] = min(abs(mapDataA.positionGrid(:,2) - desiredPos));
actualPos = mapDataA.positionGrid(index,2);

index = find(mapDataA.positionGrid(:,2) == actualPos);

profilePositionGrid = mapDataA.positionGrid(index,:);
profileDensityA = mapDataA.density.mean(index);
profileDensityB = mapDataB.density.mean(index);
profileDensityC = mapDataC.density.mean(index);

% Initialise Figure
fig = fig + 1;
set(figure(fig), 'name', 'Density_Profile', 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880])
set(gca, 'positionConstraint', 'outerPosition', ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 20, 'layer', 'top');
hold on;
    
% Plot Profiles
plot(profilePositionGrid(:,3), profileDensityA, 'color', graphColours(1), 'lineWidth', 2)
plot(profilePositionGrid(:,3), profileDensityB, 'color', graphColours(2), 'lineWidth', 2)
plot(profilePositionGrid(:,3), profileDensityC, 'color', graphColours(3), 'lineWidth', 2)
    
% Format Figure
title('.', 'color', ([254, 254, 254] / 255));
subtitle('y = - 0.1670');
axis on;
box on;
grid off;
xlim([0; 0.6]);
ylim([0; 1]);
tickData = (0.1:0.1:0.5);
xticks(tickData);
tickData = (0.2:0.2:0.8);
yticks(tickData);
xlabel('{$z_{\ell}$}', 'interpreter', 'latex');
ylabel('{Normalised Spray Density}', 'interpreter', 'latex');
xtickformat('%+.2g');
ytickformat('%+.2g');
ztickformat('%+.2g');
hold off;

tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
hold off;

pause(1);


%%

desiredPos = (0.167 / 1.044);
[~, index] = min(abs(mapDataA.positionGrid(:,2) - desiredPos));
actualPos = mapDataA.positionGrid(index,2);

index = find(mapDataA.positionGrid(:,2) == actualPos);

profilePositionGrid = mapDataA.positionGrid(index,:);
profileDensityA = mapDataA.density.mean(index);
profileDensityB = mapDataB.density.mean(index);
profileDensityC = mapDataC.density.mean(index);

% Initialise Figure
fig = fig + 1;
set(figure(fig), 'name', 'Density_Profile', 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880])
set(gca, 'positionConstraint', 'outerPosition', ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 20, 'layer', 'top');
hold on;
    
% Plot Profiles
plot(profilePositionGrid(:,3), profileDensityA, 'color', graphColours(1), 'lineWidth', 2)
plot(profilePositionGrid(:,3), profileDensityB, 'color', graphColours(2), 'lineWidth', 2)
plot(profilePositionGrid(:,3), profileDensityC, 'color', graphColours(3), 'lineWidth', 2)
    
% Format Figure
title('.', 'color', ([254, 254, 254] / 255));
subtitle('y = +0.1670');
axis on;
box on;
grid off;
xlim([0; 0.6]);
ylim([0; 1]);
tickData = (0.1:0.1:0.5);
xticks(tickData);
tickData = (0.2:0.2:0.8);
yticks(tickData);
xlabel('{$z_{\ell}$}', 'interpreter', 'latex');
ylabel('{Normalised Spray Density}', 'interpreter', 'latex');
xtickformat('%+.2g');
ytickformat('%+.2g');
ztickformat('%+.2g');
hold off;

tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
hold off;

pause(1);


%%

desiredPos = (0.1945 / 1.044);
[~, index] = min(abs(mapDataA.positionGrid(:,3) - desiredPos));
actualPos = mapDataA.positionGrid(index,3);

index = find(mapDataA.positionGrid(:,3) == actualPos);

profilePositionGrid = mapDataA.positionGrid(index,:);
profileDensityA = mapDataA.density.mean(index);
profileDensityB = mapDataB.density.mean(index);
profileDensityC = mapDataC.density.mean(index);

% Initialise Figure
fig = fig + 1;
set(figure(fig), 'name', 'Density_Profile', 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880])
set(gca, 'positionConstraint', 'outerPosition', ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 20, 'layer', 'top');
hold on;
    
% Plot Profiles
plot(profilePositionGrid(:,2), profileDensityA, 'color', graphColours(1), 'lineWidth', 2)
plot(profilePositionGrid(:,2), profileDensityB, 'color', graphColours(2), 'lineWidth', 2)
plot(profilePositionGrid(:,2), profileDensityC, 'color', graphColours(3), 'lineWidth', 2)
    
% Format Figure
title('.', 'color', ([254, 254, 254] / 255));
subtitle('z = +0.1945');
axis on;
box on;
grid off;
xlim([-0.4; 0.2]);
ylim([0; 1]);
tickData = (-0.3:0.1:0.1);
xticks(tickData);
tickData = (0.2:0.2:0.8);
yticks(tickData);
xlabel('{$y_{\ell}$}', 'interpreter', 'latex');
ylabel('{Normalised Spray Density}', 'interpreter', 'latex');
xtickformat('%+.2g');
ytickformat('%+.2g');
ztickformat('%+.2g');
hold off;

tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
hold off;

pause(1);