%% Velocity Plots v1.0

clearvars;
close all;
clc;

fig = 0;
figHold = 0; %#ok<*NASGU>

disp ('===================');
disp ('Velocity Plots v1.0');
disp ('===================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit


%% Case Selection

disp('CASE SELECTION');
disp('--------------');
disp(' ');

disp('Select Case:');
% caseFolder = uigetdir('~/Documents/Engineering/PhD/Data/Numerical/OpenFOAM');
caseFolder = ('~/Documents/Engineering/PhD/Data/Numerical/OpenFOAM/Windsor_Square_wW/kwSSTDES');
disp(['Case: ', caseFolder]);
disp(' ');
disp(' ');

if contains(caseFolder, 'Windsor_Square')
	xDims = [-0.56075; 0.48325] / 1.044;
	yDims = [-(0.389 / 2); (0.389 / 2)] / 1.044;
	zDims = [0.05; 0.339] / 1.044;

	xLims = [0.484; 1.1] / 1.044;
	yLims = [-0.207; 0.207] / 1.044;
	zLims = [0.0375; 0.3515] / 1.044;
else
    error('Unsupported Case')
end


%% x = 0.630 m (CFD)

if ~exist([caseFolder, '/velocityDataX628.csv'], 'file')
    error('No y-Plane Velocity Data Found for Target Case');
end

import = importdata([caseFolder, '/velocityDataX628.csv']);

flowFieldData.CFD.x630.points = import.data(:,4:6) / 1.044;
flowFieldData.CFD.x630.velocity = import.data(:,1:3) / 40;

[flowFieldData.CFD.x630.points, index, ~] = unique(flowFieldData.CFD.x630.points, 'rows', 'stable');
flowFieldData.CFD.x630.velocity = flowFieldData.CFD.x630.velocity(index,1:3);
[flowFieldData.CFD.x630.velocity, index, ~] = unique(flowFieldData.CFD.x630.velocity, 'rows', 'stable');
flowFieldData.CFD.x630.points = flowFieldData.CFD.x630.points(index,1:3);

[y, z] = meshgrid(min(yLims)-0.01:0.001:max(yLims)+0.01, min(zLims)-0.01:0.001:max(zLims)+0.01);

uInterp = scatteredInterpolant(flowFieldData.CFD.x630.points(:,2), flowFieldData.CFD.x630.points(:,3), flowFieldData.CFD.x630.velocity(:,1));
u = uInterp(y, z);

vInterp = scatteredInterpolant(flowFieldData.CFD.x630.points(:,2), flowFieldData.CFD.x630.points(:,3), flowFieldData.CFD.x630.velocity(:,2));
v = vInterp(y, z);

wInterp = scatteredInterpolant(flowFieldData.CFD.x630.points(:,2), flowFieldData.CFD.x630.points(:,3), flowFieldData.CFD.x630.velocity(:,3));
w = wInterp(y, z);

velMag = sqrt(u.^2 + v.^2 + w.^2);

% Figure Setup
fig = fig + 1;
figure('name', 'x = 0.630 m (CFD)');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
contourf(y, z, velMag, 16, 'edgeColor', 'none');
even_stream_line(y, z, v, w, 2, 4, 'color', 'k');

% Figure Formatting
xlim([min(yLims), max(yLims)]);
ylim([min(zLims), max(zLims)]);
tickData = min(xlim):((max(xlim) - min(xlim)) / 6):max(xlim);
xticks(tickData(2:end-1));
tickData = min(ylim):((max(ylim) - min(ylim)) / 6):max(ylim);
yticks(tickData(2:end-1));
xtickformat('%.3f');
ytickformat('%.3f');
caxis([0, 1]);
xlabel({' ', 'y (\it{l})'});
ylabel({'z (\it{l})', ' '});
box on;
colormap viridis;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
		 'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

namePos = max(strfind(caseFolder, '/'));
savefig(fig, ['~/MATLAB/Figures/', caseFolder(namePos(end):end), '_Velocity_x630']);
print(fig, ['~/MATLAB/Figures/', caseFolder(namePos(end):end), '_Velocity_x630'], '-dpng', '-r300');

clearvars -except fig figHold caseFolder xDims yDims zDims xLims yLims zLims flowFieldData;


%% x = 0.922 m (CFD)

if ~exist([caseFolder, '/velocityDataX918.csv'], 'file')
    error('No y-Plane Velocity Data Found for Target Case');
end

import = importdata([caseFolder, '/velocityDataX918.csv']);

flowFieldData.CFD.x922.points = import.data(:,4:6) / 1.044;
flowFieldData.CFD.x922.velocity = import.data(:,1:3) / 40;

[flowFieldData.CFD.x922.points, index, ~] = unique(flowFieldData.CFD.x922.points, 'rows', 'stable');
flowFieldData.CFD.x922.velocity = flowFieldData.CFD.x922.velocity(index,1:3);
[flowFieldData.CFD.x922.velocity, index, ~] = unique(flowFieldData.CFD.x922.velocity, 'rows', 'stable');
flowFieldData.CFD.x922.points = flowFieldData.CFD.x922.points(index,1:3);

[y, z] = meshgrid(min(yLims)-0.01:0.001:max(yLims)+0.01, min(zLims)-0.01:0.001:max(zLims)+0.01);

uInterp = scatteredInterpolant(flowFieldData.CFD.x922.points(:,2), flowFieldData.CFD.x922.points(:,3), flowFieldData.CFD.x922.velocity(:,1));
u = uInterp(y, z);

vInterp = scatteredInterpolant(flowFieldData.CFD.x922.points(:,2), flowFieldData.CFD.x922.points(:,3), flowFieldData.CFD.x922.velocity(:,2));
v = vInterp(y, z);

wInterp = scatteredInterpolant(flowFieldData.CFD.x922.points(:,2), flowFieldData.CFD.x922.points(:,3), flowFieldData.CFD.x922.velocity(:,3));
w = wInterp(y, z);

velMag = sqrt(u.^2 + v.^2 + w.^2);

% Figure Setup
fig = fig + 1;
figure('name', 'x = 0.922 m (CFD)');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
contourf(y, z, velMag, 16, 'edgeColor', 'none');
even_stream_line(y, z, v, w, 2, 4, 'color', 'k');

% Figure Formatting
xlim([min(yLims), max(yLims)]);
ylim([min(zLims), max(zLims)]);
tickData = min(xlim):((max(xlim) - min(xlim)) / 6):max(xlim);
xticks(tickData(2:end-1));
tickData = min(ylim):((max(ylim) - min(ylim)) / 6):max(ylim);
yticks(tickData(2:end-1));
xtickformat('%.3f');
ytickformat('%.3f');
caxis([0, 1]);
xlabel({' ', 'y (\it{l})'});
ylabel({'z (\it{l})', ' '});
box on;
colormap viridis;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
		 'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

namePos = max(strfind(caseFolder, '/'));
savefig(fig, ['~/MATLAB/Figures/', caseFolder(namePos(end):end), '_Velocity_x922']);
print(fig, ['~/MATLAB/Figures/', caseFolder(namePos(end):end), '_Velocity_x922'], '-dpng', '-r300');

clearvars -except fig figHold caseFolder xDims yDims zDims xLims yLims zLims flowFieldData;


%% y = 0 m (CFD)

if ~exist([caseFolder, '/velocityDataY0.csv'], 'file')
    error('No y-Plane Velocity Data Found for Target Case');
end

import = importdata([caseFolder, '/velocityDataY0.csv']);

flowFieldData.CFD.y0.points = import.data(:,4:6) / 1.044;
flowFieldData.CFD.y0.velocity = import.data(:,1:3) / 40;

[flowFieldData.CFD.y0.points, index, ~] = unique(flowFieldData.CFD.y0.points, 'rows', 'stable');
flowFieldData.CFD.y0.velocity = flowFieldData.CFD.y0.velocity(index,1:3);
[flowFieldData.CFD.y0.velocity, index, ~] = unique(flowFieldData.CFD.y0.velocity, 'rows', 'stable');
flowFieldData.CFD.y0.points = flowFieldData.CFD.y0.points(index,1:3);

[x, z] = meshgrid(min(xLims)-0.01:0.001:max(xLims)+0.01, min(zLims)-0.01:0.001:max(zLims)+0.01);

uInterp = scatteredInterpolant(flowFieldData.CFD.y0.points(:,1), flowFieldData.CFD.y0.points(:,3), flowFieldData.CFD.y0.velocity(:,1));
u = uInterp(x, z);

wInterp = scatteredInterpolant(flowFieldData.CFD.y0.points(:,1), flowFieldData.CFD.y0.points(:,3), flowFieldData.CFD.y0.velocity(:,3));
w = wInterp(x, z);

velMag = sqrt(u.^2 + w.^2);

% Figure Setup
fig = fig + 1;
figure('name', 'y = 0 m (CFD)');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
contourf(x, z, velMag, 16, 'edgeColor', 'none');
even_stream_line(x, z, u, w, 2, 4, 'color', 'k');
scatter(0.7055, 0.2909, 50, [0.71765, 0.00000, 0.38431]);
scatter(0.6685, 0.0899, 50, [0.71765, 0.00000, 0.38431]);
rectangle('position', ([0.405, 0.05, 0.07825, 0.289] / 1.044), 'faceColor', [0.5, 0.5, 0.5], 'edgeColor', 'none');

% Figure Formatting
xlim([(0.44 / 1.044), max(xLims)]);
ylim([min(zLims), max(zLims)]);
tickData = min(xlim):((max(xlim) - min(xlim)) / 6):max(xlim);
xticks(tickData(2:end-1));
tickData = min(ylim):((max(ylim) - min(ylim)) / 6):max(ylim);
yticks(tickData(2:end-1));
xtickformat('%.3f');
ytickformat('%.3f');
caxis([0, 1]);
xlabel({' ', 'x (\it{l})'});
ylabel({'z (\it{l})', ' '});
box on;
colormap viridis;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
		 'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

namePos = max(strfind(caseFolder, '/'));
savefig(fig, ['~/MATLAB/Figures/', caseFolder(namePos(end):end), '_Velocity_y0']);
print(fig, ['~/MATLAB/Figures/', caseFolder(namePos(end):end), '_Velocity_y0'], '-dpng', '-r300');

clearvars -except fig figHold caseFolder xDims yDims zDims xLims yLims zLims flowFieldData;


%% z = 0.194 m (CFD)

if ~exist([caseFolder, '/velocityDataZ194.csv'], 'file')
    error('No z-Plane Velocity Data Found for Target Case');
end

import = importdata([caseFolder, '/velocityDataZ194.csv']);

flowFieldData.CFD.z194.points = import.data(:,4:6) / 1.044;
flowFieldData.CFD.z194.velocity = import.data(:,1:3) / 40;

[flowFieldData.CFD.z194.points, index, ~] = unique(flowFieldData.CFD.z194.points, 'rows', 'stable');
flowFieldData.CFD.z194.velocity = flowFieldData.CFD.z194.velocity(index,1:3);
[flowFieldData.CFD.z194.velocity, index, ~] = unique(flowFieldData.CFD.z194.velocity, 'rows', 'stable');
flowFieldData.CFD.z194.points = flowFieldData.CFD.z194.points(index,1:3);

[x, y] = meshgrid(min(xLims)-0.01:0.001:max(xLims)+0.01, min(yLims)-0.01:0.001:max(yLims)+0.01);

uInterp = scatteredInterpolant(flowFieldData.CFD.z194.points(:,1), flowFieldData.CFD.z194.points(:,2), flowFieldData.CFD.z194.velocity(:,1));
u = uInterp(x, y);

vInterp = scatteredInterpolant(flowFieldData.CFD.z194.points(:,1), flowFieldData.CFD.z194.points(:,2), flowFieldData.CFD.z194.velocity(:,2));
v = vInterp(x, y);

velMag = sqrt(u.^2 + v.^2);

% Figure Setup
fig = fig + 1;
figure('name', 'z = 0.194 m (CFD)');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
contourf(x, y, velMag, 16, 'edgeColor', 'none');
even_stream_line(x, y, u, v, 2, 4, 'color', 'k');
scatter(0.5895, 0.1307, 50, [0.71765, 0.00000, 0.38431]);
scatter(0.5614, -0.1287, 50, [0.71765, 0.00000, 0.38431]);
rectangle('position', ([0.405, -0.1945, 0.07825, 0.389] / 1.044), 'faceColor', [0.5, 0.5, 0.5], 'edgeColor', 'none');

% Figure Formatting
xlim([(0.44 / 1.044), max(xLims)]);
ylim([min(yLims), max(yLims)]);
tickData = min(xlim):((max(xlim) - min(xlim)) / 6):max(xlim);
xticks(tickData(2:end-1));
tickData = min(ylim):((max(ylim) - min(ylim)) / 6):max(ylim);
yticks(tickData(2:end-1));
xtickformat('%.3f');
ytickformat('%.3f');
caxis([0, 1]);
xlabel({' ', 'x (\it{l})'});
ylabel({'y (\it{l})', ' '});
box on;
colormap viridis;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
		 'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

namePos = max(strfind(caseFolder, '/'));
savefig(fig, ['~/MATLAB/Figures/', caseFolder(namePos(end):end), '_Velocity_z194']);
print(fig, ['~/MATLAB/Figures/', caseFolder(namePos(end):end), '_Velocity_z194'], '-dpng', '-r300');

clearvars -except fig figHold caseFolder xDims yDims zDims xLims yLims zLims flowFieldData;


%% x = 0.630 m (Exp)

import = importdata('~/Documents/Engineering/PhD/Data/Experimental/Windsor Experimental Dataset/Set A - Baseline Flow/FlowField_Mean/SquareBack_WW_0_0yaw_X=0630m_FlowField.csv');
import.data(import.data == -9999) = NaN; % Prevents plotting of invalid vectors

flowFieldData.exp.x630.points = import.data(:,1:3) / 1.044;
flowFieldData.exp.x630.velocity = import.data(:,4:6) / 40;

[flowFieldData.exp.x630.points, index, ~] = unique(flowFieldData.exp.x630.points, 'rows', 'stable');
flowFieldData.exp.x630.velocity = flowFieldData.exp.x630.velocity(index,1:3);
[flowFieldData.exp.x630.velocity, index, ~] = unique(flowFieldData.exp.x630.velocity, 'rows', 'stable');
flowFieldData.exp.x630.points = flowFieldData.exp.x630.points(index,1:3);

[y, z] = meshgrid(min(yLims)-0.01:0.001:max(yLims)+0.01, min(zLims)-0.01:0.001:max(zLims)+0.01);

uInterp = scatteredInterpolant(flowFieldData.exp.x630.points(:,2), flowFieldData.exp.x630.points(:,3), flowFieldData.exp.x630.velocity(:,1));
u = uInterp(y, z);

vInterp = scatteredInterpolant(flowFieldData.exp.x630.points(:,2), flowFieldData.exp.x630.points(:,3), flowFieldData.exp.x630.velocity(:,2));
v = vInterp(y, z);

wInterp = scatteredInterpolant(flowFieldData.exp.x630.points(:,2), flowFieldData.exp.x630.points(:,3), flowFieldData.exp.x630.velocity(:,3));
w = wInterp(y, z);

velMag = sqrt(u.^2 + v.^2 + w.^2);

% Figure Setup
fig = fig + 1;
figure('name', 'x = 0.630 m (Exp)');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
contourf(y, z, velMag, 16, 'edgeColor', 'none');
even_stream_line(y, z, v, w, 2, 4, 'color', 'k');

% Figure Formatting
xlim([min(yLims), max(yLims)]);
ylim([min(zLims), max(zLims)]);
tickData = min(xlim):((max(xlim) - min(xlim)) / 6):max(xlim);
xticks(tickData(2:end-1));
tickData = min(ylim):((max(ylim) - min(ylim)) / 6):max(ylim);
yticks(tickData(2:end-1));
xtickformat('%.3f');
ytickformat('%.3f');
caxis([0, 1]);
xlabel({' ', 'y (\it{l})'});
ylabel({'z (\it{l})', ' '});
box on;
colormap viridis;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
		 'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

savefig(fig, '~/MATLAB/Figures/Windsor_Square_wW_Exp_Velocity_x630');
print(fig, '~/MATLAB/Figures/Windsor_Square_wW__Exp_Velocity_x630', '-dpng', '-r300');

clearvars -except fig figHold caseFolder xDims yDims zDims xLims yLims zLims flowFieldData;


%% x = 0.922 m (Exp)

import = importdata('~/Documents/Engineering/PhD/Data/Experimental/Windsor Experimental Dataset/Set A - Baseline Flow/FlowField_Mean/SquareBack_WW_0_0yaw_X=0922m_FlowField.csv');
import.data(import.data == -9999) = NaN; % Prevents plotting of invalid vectors

flowFieldData.exp.x922.points = import.data(:,1:3) / 1.044;
flowFieldData.exp.x922.velocity = import.data(:,4:6) / 40;

[flowFieldData.exp.x922.points, index, ~] = unique(flowFieldData.exp.x922.points, 'rows', 'stable');
flowFieldData.exp.x922.velocity = flowFieldData.exp.x922.velocity(index,1:3);
[flowFieldData.exp.x922.velocity, index, ~] = unique(flowFieldData.exp.x922.velocity, 'rows', 'stable');
flowFieldData.exp.x922.points = flowFieldData.exp.x922.points(index,1:3);

[y, z] = meshgrid(min(yLims)-0.01:0.001:max(yLims)+0.01, min(zLims)-0.01:0.001:max(zLims)+0.01);

uInterp = scatteredInterpolant(flowFieldData.exp.x922.points(:,2), flowFieldData.exp.x922.points(:,3), flowFieldData.exp.x922.velocity(:,1));
u = uInterp(y, z);

vInterp = scatteredInterpolant(flowFieldData.exp.x922.points(:,2), flowFieldData.exp.x922.points(:,3), flowFieldData.exp.x922.velocity(:,2));
v = vInterp(y, z);

wInterp = scatteredInterpolant(flowFieldData.exp.x922.points(:,2), flowFieldData.exp.x922.points(:,3), flowFieldData.exp.x922.velocity(:,3));
w = wInterp(y, z);

velMag = sqrt(u.^2 + v.^2 + w.^2);

% Figure Setup
fig = fig + 1;
figure('name', 'x = 0.922 m (Exp)');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
contourf(y, z, velMag, 16, 'edgeColor', 'none');
even_stream_line(y, z, v, w, 2, 4, 'color', 'k');

% Figure Formatting
xlim([min(yLims), max(yLims)]);
ylim([min(zLims), max(zLims)]);
tickData = min(xlim):((max(xlim) - min(xlim)) / 6):max(xlim);
xticks(tickData(2:end-1));
tickData = min(ylim):((max(ylim) - min(ylim)) / 6):max(ylim);
yticks(tickData(2:end-1));
xtickformat('%.3f');
ytickformat('%.3f');
caxis([0, 1]);
xlabel({' ', 'y (\it{l})'});
ylabel({'z (\it{l})', ' '});
box on;
colormap viridis;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
		 'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

savefig(fig, '~/MATLAB/Figures/Windsor_Square_wW_Exp_Velocity_x922');
print(fig, '~/MATLAB/Figures/Windsor_Square_wW__Exp_Velocity_x922', '-dpng', '-r300');

clearvars -except fig figHold caseFolder xDims yDims zDims xLims yLims zLims flowFieldData;


%% y = 0 m (Exp)

import = importdata('~/Documents/Engineering/PhD/Data/Experimental/Windsor Experimental Dataset/Set A - Baseline Flow/FlowField_Mean/SquareBack_WW_0_0yaw_Y=0m_FlowField.csv');
import.data(import.data == -9999) = NaN; % Prevents plotting of invalid vectors

flowFieldData.exp.y0.points = import.data(:,1:3) / 1.044;
flowFieldData.exp.y0.velocity = import.data(:,4:6) / 40;

[flowFieldData.exp.y0.points, index, ~] = unique(flowFieldData.exp.y0.points, 'rows', 'stable');
flowFieldData.exp.y0.velocity = flowFieldData.exp.y0.velocity(index,1:3);
[flowFieldData.exp.y0.velocity, index, ~] = unique(flowFieldData.exp.y0.velocity, 'rows', 'stable');
flowFieldData.exp.y0.points = flowFieldData.exp.y0.points(index,1:3);

[x, z] = meshgrid(min(xLims)-0.01:0.001:max(xLims)+0.01, min(zLims)-0.01:0.001:max(zLims)+0.01);

uInterp = scatteredInterpolant(flowFieldData.exp.y0.points(:,1), flowFieldData.exp.y0.points(:,3), flowFieldData.exp.y0.velocity(:,1));
u = uInterp(x, z);

wInterp = scatteredInterpolant(flowFieldData.exp.y0.points(:,1), flowFieldData.exp.y0.points(:,3), flowFieldData.exp.y0.velocity(:,3));
w = wInterp(x, z);

velMag = sqrt(u.^2 + w.^2);

% Figure Setup
fig = fig + 1;
figure('name', 'y = 0 m (Exp)');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
contourf(x, z, velMag, 16, 'edgeColor', 'none');
even_stream_line(x, z, u, w, 2, 4, 'color', 'k');
scatter(0.7055, 0.2909, 50, [0.71765, 0.00000, 0.38431]);
scatter(0.6685, 0.0899, 50, [0.71765, 0.00000, 0.38431]);
rectangle('position', ([0.405, 0.05, 0.07825, 0.289] / 1.044), 'faceColor', [0.5, 0.5, 0.5], 'edgeColor', 'none');

% Figure Formatting
xlim([(0.44 / 1.044), max(xLims)]);
ylim([min(zLims), max(zLims)]);
tickData = min(xlim):((max(xlim) - min(xlim)) / 6):max(xlim);
xticks(tickData(2:end-1));
tickData = min(ylim):((max(ylim) - min(ylim)) / 6):max(ylim);
yticks(tickData(2:end-1));
xtickformat('%.3f');
ytickformat('%.3f');
caxis([0, 1]);
xlabel({' ', 'x (\it{l})'});
ylabel({'z (\it{l})', ' '});
box on;
colormap viridis;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
		 'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

savefig(fig, '~/MATLAB/Figures/Windsor_Square_wW_Exp_Velocity_y0');
print(fig, '~/MATLAB/Figures/Windsor_Square_wW__Exp_Velocity_y0', '-dpng', '-r300');

clearvars -except fig figHold caseFolder xDims yDims zDims xLims yLims zLims flowFieldData;


%% z = 0.194 m (Exp)

import = importdata('~/Documents/Engineering/PhD/Data/Experimental/Windsor Experimental Dataset/Set A - Baseline Flow/FlowField_Mean/SquareBack_WW_0_0yaw_Z=0194m_FlowField.csv');
import.data(import.data == -9999) = NaN; % Prevents plotting of invalid vectors

flowFieldData.exp.z194.points = import.data(:,1:3) / 1.044;
flowFieldData.exp.z194.velocity = import.data(:,4:6) / 40;

[flowFieldData.exp.z194.points, index, ~] = unique(flowFieldData.exp.z194.points, 'rows', 'stable');
flowFieldData.exp.z194.velocity = flowFieldData.exp.z194.velocity(index,1:3);
[flowFieldData.exp.z194.velocity, index, ~] = unique(flowFieldData.exp.z194.velocity, 'rows', 'stable');
flowFieldData.exp.z194.points = flowFieldData.exp.z194.points(index,1:3);

[x, y] = meshgrid(min(xLims)-0.01:0.001:max(xLims)+0.01, min(yLims)-0.01:0.001:max(yLims)+0.01);

uInterp = scatteredInterpolant(flowFieldData.exp.z194.points(:,1), flowFieldData.exp.z194.points(:,2), flowFieldData.exp.z194.velocity(:,1));
u = uInterp(x, y);

vInterp = scatteredInterpolant(flowFieldData.exp.z194.points(:,1), flowFieldData.exp.z194.points(:,2), flowFieldData.exp.z194.velocity(:,2));
v = vInterp(x, y);

velMag = sqrt(u.^2 + v.^2);

% Figure Setup
fig = fig + 1;
figure('name', 'z = 0.194 m (Exp)');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
contourf(x, y, velMag, 16, 'edgeColor', 'none');
even_stream_line(x, y, u, v, 2, 4, 'color', 'k');
scatter(0.5895, 0.1307, 50, [0.71765, 0.00000, 0.38431]);
scatter(0.5614, -0.1287, 50, [0.71765, 0.00000, 0.38431]);
rectangle('position', ([0.405, -0.1945, 0.07825, 0.389] / 1.044), 'faceColor', [0.5, 0.5, 0.5], 'edgeColor', 'none');

% Figure Formatting
xlim([(0.44 / 1.044), max(xLims)]);
ylim([min(yLims), max(yLims)]);
tickData = min(xlim):((max(xlim) - min(xlim)) / 6):max(xlim);
xticks(tickData(2:end-1));
tickData = min(ylim):((max(ylim) - min(ylim)) / 6):max(ylim);
yticks(tickData(2:end-1));
xtickformat('%.3f');
ytickformat('%.3f');
caxis([0, 1]);
xlabel({' ', 'x (\it{l})'});
ylabel({'y (\it{l})', ' '});
box on;
colormap viridis;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
		 'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

savefig(fig, '~/MATLAB/Figures/Windsor_Square_wW_Exp_Velocity_z194');
print(fig, '~/MATLAB/Figures/Windsor_Square_wW_Exp_Velocity_z194', '-dpng', '-r300');

clearvars -except fig figHold flowFieldData;


%% Colour Bar

% Figure Setup
fig = fig + 1;
figure('name', 'Colour Bar');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Figure Formatting
caxis([0, 1]);
axis off;
colormap viridis;
c = colorbar('ticks', (0.1:0.2:0.9), 'location', 'south', 'axisLocation', 'out');
c.Label.String = {' ', '\it{U}'};
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
		 'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

print(fig, '~/MATLAB/Figures/Velocity_Bar', '-dpng', '-r300');

clearvars -except flowFieldData;