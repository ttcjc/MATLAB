%% LES Content Plots v1.0

clearvars;
close all;
clc;

fig = 0;
figHold = 0; %#ok<*NASGU>

disp ('======================');
disp ('LES Content Plots v1.0');
disp ('======================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit


%% Case Selection

disp('CASE SELECTION');
disp('--------------');
disp(' ');

disp('Select Case:');
caseFolder = uigetdir('~/Documents/Engineering/PhD/Data/Numerical/OpenFOAM');
disp(['Case: ', caseFolder]);
disp(' ');
disp(' ');

if contains(caseFolder, 'Windsor_Square')
    xDims = [-0.56075 ; 0.48325] / 1.044;
    yDims = [-(0.389 / 2) ; (0.389 / 2)] / 1.044;
    zDims = [0.05 ; 0.339] / 1.044;
else
    error('Unsupported Case')
end


%% y = 0 m (CFD)

if ~exist([caseFolder, '/LEScontentY0.csv'], 'file')
    error('No y-Plane LES Content Data Found for Target Case');
end

import = importdata([caseFolder, '/LEScontentY0.csv']);
data.LES.content = import.data(:,1) * 100;
data.points = import.data(:,2:4) / 1.044;

[data.LES.content, index, ~] = unique(data.LES.content, 'rows', 'stable');
data.points = data.points(index,1:3);
[data.points, index, ~] = unique(data.points, 'rows', 'stable');
data.LES.content = data.LES.content(index,1);

index = data.points(:,1) < max(xDims);
data.points(index,:) = [];
data.LES.content(index,:) = [];

contentInterp = scatteredInterpolant(data.points(:,1), data.points(:,3), data.LES.content(:,1));

[x, z] = meshgrid(min(data.points(:,1)):0.001:max(data.points(:,1)), min(data.points(:,3)):0.001:max(data.points(:,3)));

content = contentInterp(x, z);

% Figure Setup
fig = fig + 1;
figure('name', 'y = 0 (CFD)');
hold on;
set(figure(fig), 'outerPosition', [1945, 25, 750, 750]);

% Plot
contourf(x, z, content, 16, 'edgeColor', 'none');
rectangle('position', ([0.405, 0.05, 0.07825, 0.289] / 1.044), 'faceColor', [0.5, 0.5, 0.5], 'edgeColor', 'none');

% Figure Formatting
xlim([0.44, 1]);
ylim([0.035, 0.335]);
tickData = min(xlim):((max(xlim) - min(xlim)) / 6):max(xlim);
xticks(tickData(2:end-1));
tickData = min(zDims):((max(zDims) - min(zDims)) / 6):max(zDims);
yticks(tickData(2:end-1));
xtickformat('%.3f');
ytickformat('%.3f');
caxis([0, 100]);
xlabel({' ', 'x (\it{l})'});
ylabel({'z (\it{l})', ' '});
box on;
colormap viridis;
set(gca, 'units', 'normalized', 'position', [0.125, 0.125, 0.75, 0.75], ...
         'fontName', 'LM Roman 12', 'fontSize', 11, 'layer', 'top');
hold off;

print(fig, '~/MATLAB/Output/Figures/LES_Content_y0', '-dpng', '-r300');

clearvars -except fig figHold caseFolder xDims yDims zDims;


%% y = 0.194 m (CFD)

if ~exist([caseFolder, '/LEScontentZ194.csv'], 'file')
    error('No y-Plane LES Content Data Found for Target Case');
end

import = importdata([caseFolder, '/LEScontentZ194.csv']);
data.LES.content = import.data(:,1) * 100;
data.points = import.data(:,2:4) / 1.044;

[data.LES.content, index, ~] = unique(data.LES.content, 'rows', 'stable');
data.points = data.points(index,1:3);
[data.points, index, ~] = unique(data.points, 'rows', 'stable');
data.LES.content = data.LES.content(index,1);

index = data.points(:,1) < max(xDims);
data.points(index,:) = [];
data.LES.content(index,:) = [];

contentInterp = scatteredInterpolant(data.points(:,1), data.points(:,2), data.LES.content(:,1));

[x, y] = meshgrid(min(data.points(:,1)):0.001:max(data.points(:,1)), min(data.points(:,2)):0.001:max(data.points(:,2)));

content = contentInterp(x, y);

% Figure Setup
fig = fig + 1;
figure('name', 'z = 0.194 (CFD)');
hold on;
set(figure(fig), 'outerPosition', [1945, 100, 750, 750]);

% Plot
contourf(x, y, content, 16, 'edgeColor', 'none');
rectangle('position', ([0.405, -0.1945, 0.07825, 0.389] / 1.044), 'faceColor', [0.5, 0.5, 0.5], 'edgeColor', 'none');

% Figure Formatting
xlim([0.44, 1]);
ylim([-0.2, 0.2]);
tickData = min(xlim):((max(xlim) - min(xlim)) / 6):max(xlim);
xticks(tickData(2:end-1));
tickData = min(yDims):((max(yDims) - min(yDims)) / 6):max(yDims);
yticks(tickData(2:end-1));
xtickformat('%.3f');
ytickformat('%.3f');
caxis([0, 100]);
xlabel({' ', 'x (\it{l})'});
ylabel({'y (\it{l})', ' '});
box on;
colormap viridis;
set(gca, 'units', 'normalized', 'position', [0.125, 0.125, 0.75, 0.75], ...
         'fontName', 'LM Roman 12', 'fontSize', 11, 'layer', 'top');

print(fig, '~/MATLAB/Output/Figures/LES_Content_z194', '-dpng', '-r300');

clearvars -except fig figHold xDims yDims zDims;


%% Colour Bar

% Figure Setup
fig = fig + 1;
figure('name', 'Colour Bar');
hold on;
set(figure(fig), 'outerPosition', [2332.5, 25, 750, 750]);

% Figure Formatting
caxis([0, 100]);
axis off;
colormap viridis;
c = colorbar('ticks', (10:20:90), 'location', 'southOutside');
c.Label.String = {' ', 'Resolved Turbulence Kinetic Energy (%)'};
set(gca, 'fontName', 'LM Roman 12', 'fontSize', 11, 'layer', 'top');
hold off;

print(fig, '~/MATLAB/Output/Figures/LES_Content_Bar', '-dpng', '-r300');


%% Cleaning

clearvars;
disp(' ');