%% Aerodynamic Coefficient Calculator v2.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

normalise = false; % Normalisation of Dimensions

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp('=======================================');
disp('Aerodynamic Coefficient Calculator v2.0');
disp('=======================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats


%% Initialise Case

% Select Case
disp('Case Selection');
disp('---------------');

caseFolder = uigetdir('~/OpenFOAM/', 'Select Case');

namePos = max(strfind(caseFolder, '/')) + 1;
caseName = caseFolder(namePos:end);

disp(' ');

disp(['Case: ', caseName]);

% Confirm Support
if ~contains(caseName, ["Run_Test", "Windsor"])
    error('Invalid Case Directory (Unsupported Case Type)');
end

disp(' ');
disp(' ');


%% Acquire Coefficient Data

% Confirm Data Availability
if exist([caseFolder, '/postProcessing/forceCoeffs'], 'dir')
    coeffDirs = dir([caseFolder, '/postProcessing/forceCoeffs']);
else
    error('Invalid Case Directory (No Coefficient Data Found)');
end

% Identify Coefficient Directories
i = 1;
while i <= height(coeffDirs)

    if isnan(str2double(coeffDirs(i).name))
        coeffDirs(i) = [];
    else
        i = i + 1;
    end

end
clear i;

% Collate Coefficient Data
coeffData.time = [];
coeffData.Cd = [];
coeffData.Cl = [];

for i = 1:height(coeffDirs)
    fileID = fopen([caseFolder, '/postProcessing/forceCoeffs/', coeffDirs(i).name, '/forceCoeffsExtended.dat']);
    content = textscan(fileID, '%f %f %f %f %f %f', 'headerLines', 9, 'delimiter', '\n', 'collectOutput', true);

    coeffData.time = [coeffData.time; content{1}(:,1)];
    coeffData.Cd = [coeffData.Cd; content{1}(:,3)];
    coeffData.Cl = [coeffData.Cl; content{1}(:,4)];

    fclose(fileID);
end

% Load Experimental Values
if strcmp(caseName, 'Windsor_SB_wW_Upstream_SC')
    coeffData.Cd_Exp = 0.350;
    coeffData.Cl_Exp = 0.134;
elseif strcmp(caseName, 'Windsor_ST_20D_wW_Upstream_SC')
    coeffData.Cd_Exp = 0.328;
    coeffData.Cl_Exp = 0.107;
elseif strcmp(caseName, 'Windsor_RSST_16D_U50_wW_Upstream_SC')
    coeffData.Cd_Exp = 0.333;
    coeffData.Cl_Exp = 0.195;
else
    coeffData.Cd_Exp = nan;
    coeffData.Cl_Exp = nan;
end

% Perform Blockage Correction
if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
    U = 40; % m/s
    rho = 1.269; % kg/m^3

    Am = (0.289 * 0.389) + (2 * (0.046 * 0.055));
    At = (2 * (0.9519083 + (3.283 * tan(atan(0.0262223 / 9.44)))) * 1.32);

    Fd = coeffData.Cd * 0.5 * rho * U^2 * Am;
    Fl = coeffData.Cl * 0.5 * rho * U^2 * Am;
    
    Ucorr = U * (At / (At - Am));
    
    coeffData.Cd = Fd / (0.5 * rho * Ucorr^2 * Am);
    coeffData.Cl = Fl / (0.5 * rho * Ucorr^2 * Am);
end


%% Calculate Time-Averaged Values

averageStart = find(coeffData.time == 1);

coeffData.Cd_Mean = round(mean(coeffData.Cd(averageStart:end)), 3);
coeffData.Cl_Mean = round(mean(coeffData.Cl(averageStart:end)), 3);

coeffData.Cd_Error = ((coeffData.Cd_Exp - coeffData.Cd_Mean) / coeffData.Cd_Exp) * 100;
coeffData.Cl_Error = ((coeffData.Cl_Exp - coeffData.Cl_Mean) / coeffData.Cl_Exp) * 100;

disp(['Time-Averaged Numerical Drag Coefficient = ', num2str(coeffData.Cd_Mean)]);
disp(['Time-Averaged Numerical Lift Coefficient = ', num2str(coeffData.Cl_Mean)]);

disp(' ');

disp(['Time-Averaged Experimental Drag Coefficient = ', num2str(coeffData.Cd_Exp)]);
disp(['Time-Averaged Experimental Lift Coefficient = ', num2str(coeffData.Cl_Exp)]);

disp(' ');

disp(['Drag Coefficient Error = ', num2str(coeffData.Cd_Error), '%']);
disp(['Lift Coefficient Error = ', num2str(coeffData.Cl_Error), '%']);


%% Plot Drag Coefficient Data

% Figure Setup
fig = fig + 1;
figName = [caseName, '_Drag_Coefficient'];
set(figure(fig), 'outerPosition', [25, 25, 1275, 850], 'name', figName);
set(gca, 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
         'fontSize', 20, 'layer', 'top');
hold on;

% Plot
plot(coeffData.time(10:end), coeffData.Cd(10:end), 'lineWidth', 1.5, 'color', ([74, 24, 99] / 255));
xline(1, 'k', 'Start of Averaging', 'lineWidth', 1.5);

% Figure Formatting
axis on;
box on;
grid off;
xlim([0; 5]);
ylim([0.25; 0.45]);
tickData = (0:0.5:5);
xticks(tickData(2:(end - 1)));
tickData = (0.25:0.025:0.45);
yticks(tickData(2:(end - 1)));
xtickformat('%.1f');
ytickformat('%+.3f');
xlabel({' ', '{\bf{Time (\it{s})}}'}, 'fontName', 'LM Roman 12');
ylabel({'{\bf{C_D}}', ' '}, 'fontName', 'LM Roman 12');
set(gca, 'outerPosition', [0.05, 0.05, 0.9, 0.9]);
hold off;

pause(2);
exportgraphics(gcf, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);


%% Plot Lift Coefficient Data

% Figure Setup
fig = fig + 1;
figName = [caseName, '_Lift_Coefficient'];
set(figure(fig), 'outerPosition', [25, 25, 1275, 850], 'name', figName);
set(gca, 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
         'fontSize', 20, 'layer', 'top');
hold on;

% Plot
plot(coeffData.time(10:end), coeffData.Cl(10:end), 'lineWidth', 1.5, 'color', ([230, 0, 126] / 255));
xline(1, 'k', 'Start of Averaging', 'lineWidth', 1.5);

% Figure Formatting
axis on;
box on;
grid off;
xlim([0; 5]);
ylim([-0.05; 0.15]);
tickData = (0:0.5:5);
xticks(tickData(2:(end - 1)));
tickData = (-0.05:0.025:0.15);
yticks(tickData(2:(end - 1)));
xtickformat('%.1f');
ytickformat('%+.3f');
xlabel({' ', '{\bf{Time (\it{s})}}'}, 'fontName', 'LM Roman 12');
ylabel({'{\bf{C_L}}', ' '}, 'fontName', 'LM Roman 12');
set(gca, 'outerPosition', [0.05, 0.05, 0.9, 0.9]);
hold off;

pause(2);
exportgraphics(gcf, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);