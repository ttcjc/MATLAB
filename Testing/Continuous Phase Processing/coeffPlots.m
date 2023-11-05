%% Aerodynamic Coefficient Plotter v2

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('=======================================');
disp ('%% Aerodynamic Coefficient Plotter v2');
disp ('=======================================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Concatenated 'coeffData' to Improve Readability
% v2.0 - Updated to Support 'forceCoeffsExtended' Functionality


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
    while i <= height(coeffDirs)

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
for i = 1:height(coeffDirs)
    coeffData.directory{i,1} = coeffDirs(i,1).name;

    fileID = fopen([caseFolder, '/postProcessing/forceCoeffs/', coeffData.directory{i,1}, '/forceCoeffsExtended.dat']);
    content = textscan(fileID, '%f %f %f %f %f %f %f', 'headerLines', 12, 'delimiter', ' ', 'multipleDelimsAsOne', 1);

    coeffData.time{i,1} = content{1,1};
    coeffData.Cl{i,1} = content{1,2};
    coeffData.Cd{i,1} = content{1,3};
    coeffData.Cs{i,1} = content{1,4};
    coeffData.Cm_p{i,1} = content{1,5};
    coeffData.Cm_y{i,1} = content{1,6};
    coeffData.Cm_r{i,1} = content{1,7};

    fclose(fileID);
end

% Concatenate Coefficient Data
time = coeffData.time{1,1};
Cl = coeffData.Cl{1,1};
Cd = coeffData.Cd{1,1};
Cs = coeffData.Cs{1,1};
Cm_p = coeffData.Cm_p{1,1};
Cm_y = coeffData.Cm_y{1,1};
Cm_r = coeffData.Cm_r{1,1};

if height(coeffData.directory) > 1

    for i = 2:height(coeffData.directory)
        time = vertcat(time, coeffData.time{i,1}); %#ok<AGROW>
        Cl = vertcat(Cl, coeffData.Cl{i,1}); %#ok<AGROW>
        Cd = vertcat(Cd, coeffData.Cd{i,1}); %#ok<AGROW>
        Cs = vertcat(Cs, coeffData.Cs{i,1}); %#ok<AGROW>
        Cm_p = vertcat(Cm_p, coeffData.Cm_p{i,1}); %#ok<AGROW>
        Cm_y = vertcat(Cm_y, coeffData.Cm_y{i,1}); %#ok<AGROW>
        Cm_r = vertcat(Cm_r, coeffData.Cm_r{i,1}); %#ok<AGROW>
    end

end

coeffData.time = time;
coeffData.Cl = Cl;
coeffData.Cd = Cd;
coeffData.Cs = Cs;
coeffData.Cm_p = Cm_p;
coeffData.Cm_y = Cm_y;
coeffData.Cm_r = Cm_r;

clear time Cl Cd Cs Cm_p Cm_y Cm_r;


%% Process Stuff

Cd_Mean = zeros(height(coeffData.time),1);

startTime = 0.25;
startSample = find(coeffData.time == startTime);

for i = 1:height(coeffData.time)
    
    if i < startSample
        Cd_Mean(i) = nan;
    else
        Cd_Mean(i) = mean(coeffData.Cd(2500:i));
    end
    
end


%% Plot

% Figure Setup
fig = fig + 1;
set(figure(fig), 'outerPosition', [25, 25, 1275, 850], 'name', 'Lift Coefficient');
set(gca, 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
         'fontSize', 20, 'layer', 'top');
hold on;

% Plot
plot(coeffData.time, coeffData.Cl, 'lineWidth', 2, 'color', ([74, 24, 99] / 255));
% plot(coeffData.time, Cl_Mean, 'lineWidth', 2, 'color', ([230, 0, 126] / 255));

% Figure Formatting
xlim('padded');
ylim([0; 0.3]);
xticks([]);
yticks([]);
xlabel({' ', 'Time'});
ylabel({'Lift Coefficient', ' '});
grid on;
box on;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
         'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;


% Figure Setup
fig = fig + 1;
set(figure(fig), 'outerPosition', [25, 25, 1275, 850], 'name', 'Drag Coefficient');
set(gca, 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
         'fontSize', 20, 'layer', 'top');
hold on;

% Plot
plot(coeffData.time, coeffData.Cd, 'lineWidth', 2, 'color', ([74, 24, 99] / 255));
plot(coeffData.time, Cd_Mean, 'lineWidth', 2, 'color', ([230, 0, 126] / 255));

% Figure Formatting
xlim('padded');
ylim([0.3,0.45]);
xticks([]);
yticks([]);
xlabel({' ', 'Time'});
ylabel({'Drag Coefficient', ' '});
grid on;
box on;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
         'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;


% Figure Setup
fig = fig + 1;
set(figure(fig), 'outerPosition', [25, 25, 1275, 850], 'name', 'Side Coefficient');
set(gca, 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
         'fontSize', 20, 'layer', 'top');
hold on;

% Plot
plot(coeffData.time, coeffData.Cs, 'lineWidth', 2, 'color', ([74, 24, 99] / 255));
yline(0, 'lineStyle', '--', 'lineWidth', 2, 'color', ([230, 0, 126] / 255));

% Figure Formatting
xlim('padded');
ylim([-0.1,0.1]);
xticks([]);
yticks([]);
xlabel({' ', 'Time'});
ylabel({'Side Coefficient', ' '});
grid on;
box on;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
         'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;