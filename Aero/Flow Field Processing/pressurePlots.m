%% Pressure Plots v1.0

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('===================');
disp ('Pressure Plots v1.0');
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
caseFolder = ('~/Documents/Engineering/PhD/Data/Numerical/OpenFOAM/Windsor_Square_wW/SADDES');
disp(['Case: ', caseFolder]);
disp(' ');
disp(' ');

if contains(caseFolder, 'Windsor_Square')
    xDims = [-0.56075 ; 0.48325] / 1.044;
    yDims = [-(0.389 / 2) ; (0.389 / 2)] / 1.044;
    zDims = [0.05 ; 0.339] / 1.044;
	
	A_Model = (0.289 * 0.389) + 2 * (0.05 * 0.055);
	A_Tunnel = (2 * (0.96 + (1.695 * tan(atan(0.01 / 3.6))))) - (4 * 0.01125); % At balance point of resolution
	A_Domain = (2 * (0.96 + (4.704 * tan(atan(0.01 / 3.6))))); % At midpoint
else
    error('Unsupported Case')
end

E_Exp = A_Model / A_Tunnel;
E_CFD = A_Model / A_Domain;


%% Base (Exp)

import = importdata('~/Documents/Engineering/PhD/Data/Experimental/Windsor Experimental Dataset/Set A - Baseline Flow/Pressure_Mean/SquareBack_Pressure_Tapping_Map.csv');
data.exp.base.points(:,1:3) = import.data(51:106,1:3);
data.exp.base.points(:,1) = data.exp.base.points(:,1) + 0.65;
data.exp.base.points(:,1:3) = data.exp.base.points(:,1:3) / (1e3 * 1.044);

import = importdata('~/Documents/Engineering/PhD/Data/Experimental/Windsor Experimental Dataset/Set A - Baseline Flow/Pressure_Mean/SquareBack_WW_0_0yaw_Pressures_Averages.csv');
data.exp.base.CpMean(:,1) = import.data(1,51:106);
data.exp.base.CpRMS(:,1) = import.data(2,51:106);

[y, z] = meshgrid(min(yDims):(0.001 / 1.044):max(yDims), min(zDims):(0.001 / 1.044):max(zDims));

CpMeanInterp = scatteredInterpolant(data.exp.base.points(:,2), data.exp.base.points(:,3), data.exp.base.CpMean(:,1));
CpMean = CpMeanInterp(y, z);

CpRMSInterp = scatteredInterpolant(data.exp.base.points(:,2), data.exp.base.points(:,3), data.exp.base.CpRMS(:,1));
CpRMS = CpRMSInterp(y, z);

CpMeanCorr = (CpMean + (2 * E_Exp)) / (1 + (2 * E_Exp));

% Figure Setup
fig = fig + 1;
figure('name', 'Base Mean Cp (Exp)');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
contourf(y, z, CpMeanCorr, 16, 'edgeColor', 'none');
scatter(data.exp.base.points(:,2), data.exp.base.points(:,3), 50, [0.71765, 0.00000, 0.38431]);

% Figure Formatting
tickData = min(yDims):((max(yDims) - min(yDims)) / 6):max(yDims);
xticks(tickData(2:end-1));
tickData = min(zDims):((max(zDims) - min(zDims)) / 6):max(zDims);
yticks(tickData(2:end-1));
xtickformat('%.3f');
ytickformat('%.3f');
caxis([-0.2007, -0.0969]);
xlabel({' ', 'y (\it{l})'});
ylabel({'z (\it{l})', ' '});
box on;
colormap viridis;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
	     'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

savefig(fig, '~/MATLAB/Figures/Windsor_Square_wW_Exp_Cp_Base');
print(fig, '~/MATLAB/Figures/Windsor_Square_wW_Exp_Cp_Base', '-dpng', '-r300');

% Figure Setup
fig = fig + 1;
figure('name', 'Base RMS Cp (Exp)');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
contourf(y, z, CpRMS, 16, 'edgeColor', 'none');
scatter(data.exp.base.points(:,2), data.exp.base.points(:,3), 50, [0.71765, 0.00000, 0.38431]);

% Figure Formatting
tickData = min(yDims):((max(yDims) - min(yDims)) / 6):max(yDims);
xticks(tickData(2:end-1));
tickData = min(zDims):((max(zDims) - min(zDims)) / 6):max(zDims);
yticks(tickData(2:end-1));
xtickformat('%.3f');
ytickformat('%.3f');
% caxis([-0.2007, -0.0969]);
xlabel({' ', 'y (\it{l})'});
ylabel({'z (\it{l})', ' '});
box on;
colormap viridis;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
	     'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

savefig(fig, '~/MATLAB/Figures/Windsor_Square_wW_Exp_CpRMS_Base');
print(fig, '~/MATLAB/Figures/Windsor_Square_wW_Exp_CpRMS_Base', '-dpng', '-r300');

clearvars -except fig figHold caseFolder xDims yDims zDims E_Exp E_CFD data;


%% Centreline (Exp)

import = importdata('~/Documents/Engineering/PhD/Data/Experimental/Windsor Experimental Dataset/Set A - Baseline Flow/Pressure_Mean/SquareBack_Pressure_Tapping_Map.csv');
data.exp.centreline.points(:,1:3) = import.data(1:18,1:3);
data.exp.centreline.points(:,1) = data.exp.centreline.points(:,1) + 0.65;
data.exp.centreline.points(:,1:3) = data.exp.centreline.points(:,1:3) / (1e3 * 1.044);

import = importdata('~/Documents/Engineering/PhD/Data/Experimental/Windsor Experimental Dataset/Set A - Baseline Flow/Pressure_Mean/SquareBack_WW_0_0yaw_Pressures_Averages.csv');
data.exp.centreline.CpMean(:,1) = import.data(1,1:18);
data.exp.centreline.CpRMS(:,1) = import.data(2,1:18);

x = data.exp.centreline.points(:,1);
CpMean = data.exp.centreline.CpMean(:,1);

CpMeanCorr = (CpMean + (2 * E_Exp)) / (1 + (2 * E_Exp));

% Figure Setup
fig = fig + 1;
figure('name', 'Centreline Mean Cp');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
plot(x, CpMeanCorr, 'color', [0.71765, 0.00000, 0.38431]);

figHold = fig;
hold off;

clearvars -except fig figHold caseFolder xDims yDims zDims E_CFD data;


%% Base (CFD)

if ~exist([caseFolder, '/pressureDataBase.csv'], 'file')
    error('No Base Pressure Data Found for Target Case');
end

import = importdata([caseFolder, '/pressureDataBase.csv']);
data.CFD.base.pressure = import.data(:,1);
data.CFD.base.points = import.data(:,2:4) / 1.044;

[data.CFD.base.pressure, index, ~] = unique(data.CFD.base.pressure, 'rows', 'stable');
data.CFD.base.points = data.CFD.base.points(index,1:3);
[data.CFD.base.points, index, ~] = unique(data.CFD.base.points, 'rows', 'stable');
data.CFD.base.pressure = data.CFD.base.pressure(index,:);

index = data.CFD.base.points(:,1) ~= (0.48325 / 1.044);
data.CFD.base.points(index,:) = [];
data.CFD.base.pressure(index,:) = [];

data.CFD.base.pressure = ((data.CFD.base.pressure * 1.2047) - 0) / (0.5 * 1.2047 * 40^2);

CpMeanInterp = scatteredInterpolant(data.CFD.base.points(:,2), data.CFD.base.points(:,3), data.CFD.base.pressure(:,1));
[y, z] = meshgrid(min(yDims):(0.001 / 1.044):max(yDims), min(zDims):(0.001 / 1.044):max(zDims));
CpMean = CpMeanInterp(y, z);

CpMeanCorr = (CpMean + (2 * E_CFD)) / (1 + (2 * E_CFD));

% Figure Setup
fig = fig + 1;
figure('name', 'Base Mean Cp (CFD)');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
contourf(y, z, CpMeanCorr, 16, 'edgeColor', 'none');
scatter(data.exp.base.points(:,2), data.exp.base.points(:,3), 50, [0.71765, 0.00000, 0.38431]);

% Figure Formatting
tickData = min(yDims):((max(yDims) - min(yDims)) / 6):max(yDims);
xticks(tickData(2:end-1));
tickData = min(zDims):((max(zDims) - min(zDims)) / 6):max(zDims);
yticks(tickData(2:end-1));
xtickformat('%.3f');
ytickformat('%.3f');
% caxis([-0.2007, -0.0969]);
xlabel({' ', 'y (\it{l})'});
ylabel({'z (\it{l})', ' '});
box on;
colormap viridis;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
	     'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

namePos = max(strfind(caseFolder, '/'));
savefig(fig, ['~/MATLAB/Figures/', caseFolder(namePos(end):end), '_Cp_Base']);
print(fig, ['~/MATLAB/Figures/', caseFolder(namePos(end):end), '_Cp_Base'], '-dpng', '-r300');

clearvars -except fig figHold caseFolder xDims yDims zDims E_CFD data;


%% Centreline (CFD)

if ~exist([caseFolder, '/pressureDataCentreline.csv'], 'file')
    error('No Centreline Pressure Data Found for Target Case');
end

import = importdata([caseFolder, '/pressureDataCentreline.csv']);
data.CFD.centreline.pressure = import.data(:,1);
data.CFD.centreline.points = import.data(:,2:4) / 1.044;

[data.CFD.centreline.pressure, index, ~] = unique(data.CFD.centreline.pressure, 'rows', 'stable');
data.CFD.centreline.points = data.CFD.centreline.points(index,1:3);
[data.CFD.centreline.points, index, ~] = unique(data.CFD.centreline.points, 'rows', 'stable');
data.CFD.centreline.pressure = data.CFD.centreline.pressure(index,:);

index = data.CFD.centreline.points(:,3) < (0.12492 / 1.044);
data.CFD.centreline.points(index,:) = [];
data.CFD.centreline.pressure(index,:) = [];
index = data.CFD.centreline.points(:,1) > (0.42625 / 1.044);
data.CFD.centreline.points(index,:) = [];
data.CFD.centreline.pressure(index,:) = [];

[data.CFD.centreline.points, index] = sortrows(data.CFD.centreline.points,1);
data.CFD.centreline.pressure = data.CFD.centreline.pressure(index,:);

data.CFD.centreline.pressure = ((data.CFD.centreline.pressure * 1.2047) - 0) / (0.5 * 1.2047 * 40^2);

x = data.CFD.centreline.points(:,1);
CpMean = data.CFD.centreline.pressure(:,1);

CpMeanCorr = (CpMean + (2 * E_CFD)) / (1 + (2 * E_CFD));

% Figure Setup
figure(figHold);
hold on;

% Plot
plot(x, CpMeanCorr, 'color', [0.21176, 0.06667, 0.38824]);

% Figure Formatting
xlim([(min(xDims) - 0.025), (max(xDims) + 0.025)]);
ylim([-1, 1]);
tickData = min(xlim):((max(xlim) - min(xlim)) / 6):max(xlim);
xticks(tickData(2:end-1));
yticks(-0.9:0.45:0.9);
xtickformat('%.3f');
ytickformat('%.3f');
xlabel({' ', 'x (\it{l})'});
ylabel({'C_p', ' '});
grid on;
box on;
legend('Experiment (Varney)', 'Numerical', 'location', 'northOutside', 'orientation', 'vertical');
legend boxoff;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
	     'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');

namePos = max(strfind(caseFolder, '/'));
savefig(figHold, ['~/MATLAB/Figures/', caseFolder(namePos(end):end), '_Cp_Centreline']);
print(figHold, ['~/MATLAB/Figures/', caseFolder(namePos(end):end), '_Cp_Centreline'], '-dpng', '-r300');

clearvars -except fig figHold data;


%% Colour Bar

% Figure Setup
fig = fig + 1;
figure('name', 'Colour Bar');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Figure Formatting
caxis([-0.235, -0.0055]);
tickData = min(caxis):((max(caxis) - min(caxis)) / 6):max(caxis);
axis off;
colormap viridis;
c = colorbar('ticks', round(tickData(2:end-1),3), 'location', 'south', 'axisLocation', 'out');
c.Label.String = {' ', 'C_p'};
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
	     'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

print(fig, '~/MATLAB/Figures/Pressure_Bar', '-dpng', '-r300');

clearvars -except data;