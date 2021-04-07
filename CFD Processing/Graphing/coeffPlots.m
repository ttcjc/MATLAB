%% Aerodynamic Coefficient Plotter v1.1

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('=======================================');
disp ('%% Aerodynamic Coefficient Plotter v1.1');
disp ('=======================================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Concatenated 'coeffData' to Improve Readability


%% Case Selection

disp('CASE SELECTION');
disp('--------------');
disp(' ');

disp('Select Case Folder:');
caseFolder = uigetdir('~/OpenFOAM');

disp(['Case: ', caseFolder]);
disp(' ');


%% Data Acquisition

% Confirm Case Validity and Identify Time Directories
if exist([caseFolder, '/postProcessing/forceCoeffs'], 'dir')
    coeffDirs = dir([caseFolder, '/postProcessing/forceCoeffs']);

    i = 1;
    while i <= size(coeffDirs,1)

        if isnan(str2double(coeffDirs(i,1).name))
            coeffDirs(i,:) = [];
        else
            i = i + 1;
        end

    end

    if isempty(coeffDirs)
        error('No Coefficient Data Found for Target Case');
    end

else
    error('No Coefficient Data Found for Target Case');
end

% Read Coefficient Data
for i = 1:size(coeffDirs,1)
    coeffData.directory{i,1} = coeffDirs(i,1).name;

    fileID = fopen([caseFolder, '/postProcessing/forceCoeffs/', coeffData.directory{i,1}, '/forceCoeffs.dat']);
    content = textscan(fileID, '%f %f %f %f %f %f', 'headerLines', 9, 'delimiter', ' ', 'multipleDelimsAsOne', 1);

    coeffData.time{i,1} = content{1,1};
%     coeffData.Cm{i,1} = content{1,2};
    coeffData.Cd{i,1} = content{1,3};
    coeffData.Cl{i,1} = content{1,4};
%     coeffData.Cl_f{i,1} = content{1,5};
%     coeffData.Cl_r{i,1} = content{1,6};

    fclose(fileID);
end

% Concatenate Coefficient Data
time = coeffData.time{1,1};
Cd = coeffData.Cd{1,1};
Cl = coeffData.Cl{1,1};

if size(coeffData.directory,1) > 1

    for i = 2:size(coeffData.directory,1)
        time = vertcat(time, coeffData.time{i,1}); %#ok<AGROW>
        Cd = vertcat(Cd, coeffData.Cd{i,1}); %#ok<AGROW>
        Cl = vertcat(Cl, coeffData.Cl{i,1}); %#ok<AGROW>
    end

end

coeffData.time = time;
coeffData.Cd = Cd;
coeffData.Cl = Cl;

% Load experimental values
if contains(caseFolder, 'Windsor_Square_wW')
    coeffData.Cd_Exp_Mean = 0.3889;
    coeffData.Cl_Exp_Mean = 0.1435;

    U = 39.56;
    rho = 1.269;
else
    coeffData.Cd_Exp_Mean = nan;
    coeffData.Cl_Exp_Mean = nan;
    
    U = nan;
    rho = nan;
end


%% Blockage Correction

% Experimental
if ~isnan(coeffData.Cd_Exp_Mean)
    At = ((2 * (0.96 + (1.695 * tan(atan(0.01 / 3.6))))) * 1.32) - (4 * 0.01125);
    Am = (0.289 * 0.389) + (2 * (0.05 * 0.055));

    Ucorr = (U * At) / (At - Am);

    Fd = coeffData.Cd_Exp_Mean * (0.5 * rho * U^2 * Am);
    Fl = coeffData.Cl_Exp_Mean * (0.5 * rho * U^2 * Am);

    coeffData.Cd_Exp_Mean_Corr = Fd / (0.5 * rho * Ucorr^2 * Am);
    coeffData.Cl_Exp_Mean_Corr = Fl / (0.5 * rho * Ucorr^2 * Am);
else
    coeffData.Cd_Exp_Mean_Corr = nan;
    coeffData.Cl_Exp_Mean_Corr = nan;
end

% Numerical
if contains(caseFolder, 'Upstream')    
    At = (2 * (0.96 + (3.379 * tan(atan(0.02613 / 9.408)))) * 1.32);
else
    At = (2 * (0.96 + (4.704 * tan(atan(0.02613 / 9.408)))) * 1.32);
end

Am = (0.289 * 0.389) + (2 * (0.046 * 0.055));

U = 40;
rho = 1.269;

Ucorr = (U * At) / (At - Am);

Fd = coeffData.Cd * (0.5 * rho * U^2 * Am);
Fl = coeffData.Cl * (0.5 * rho * U^2 * Am);

coeffData.Cd_Corr = Fd / (0.5 * rho * Ucorr^2 * Am);
coeffData.Cl_Corr = Fl / (0.5 * rho * Ucorr^2 * Am);


%% Time-Averaged Values

coeffData.Cd_Corr_Mean = round(mean(coeffData.Cd_Corr(round(end / 3):end)),4);
coeffData.Cl_Corr_Mean = round(mean(coeffData.Cl_Corr(round(end / 3):end)),4);

coeffData.Cd_Error = ((coeffData.Cd_Exp_Mean_Corr - coeffData.Cd_Corr_Mean) / coeffData.Cd_Exp_Mean_Corr) * 100;
coeffData.Cl_Error = ((coeffData.Cl_Exp_Mean_Corr - coeffData.Cl_Corr_Mean) / coeffData.Cl_Exp_Mean_Corr) * 100;

disp(['Time-Averaged Numerical Drag Coefficient = ', num2str(coeffData.Cd_Corr_Mean)]);
disp(['Time-Averaged Numerical Lift Coefficient = ', num2str(coeffData.Cl_Corr_Mean)]);
disp(' ');
disp(['Experimental Drag Coefficient = ', num2str(coeffData.Cd_Exp_Mean_Corr)]);
disp(['Experimental Lift Coefficient = ', num2str(coeffData.Cl_Exp_Mean_Corr)]);
disp(' ');
disp(['Drag Coefficient Error = ', num2str(coeffData.Cd_Error), '%']);
disp(['Lift Coefficient Error = ', num2str(coeffData.Cl_Error), '%']);


%% Plot

% Figure Setup
fig = fig + 1;
figure('name', 'Force Coefficients');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
plot(coeffData.time, coeffData.Cd_Corr, 'color', [0.21176, 0.06667, 0.38824]);
plot(coeffData.time, coeffData.Cl_Corr, 'color', [0.71765, 0.00000, 0.38431]);

% Figure Formatting
xlim([0,4])
ylim([-0.2,0.6]);
xticks(0.8:0.8:3.2);
yticks(-0.04:0.16:0.44);
xtickformat('%.1f');
ytickformat('%.2f');
xlabel({' ', 'Time (\it{s})'});
ylabel({'Force Coefficient', ' '});
grid on;
box on;
legend('C_D', 'C_L', 'location', 'northOutside', 'orientation', 'vertical');
legend boxoff;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
         'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

namePos = max(strfind(caseFolder, '/'));
savefig(fig, ['~/MATLAB/Output/Figures/', caseFolder(namePos(end):end), '_Force_Coeffs_']);
print(fig, ['~/MATLAB/Output/Figures/', caseFolder(namePos(end):end), '_Force_Coeffs_'], '-dpng', '-r300');


%% Cleaning

clearvars -except coeffData;
disp(' ');