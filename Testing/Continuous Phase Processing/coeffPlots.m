%% Aerodynamic Coefficient Plotter v2

run preamble;

%#ok<*UNRCH>

blockageCorrection = true; % Perform Blockage Correction

disp ('=======================================');
disp ('%% Aerodynamic Coefficient Plotter v2');
disp ('=======================================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Concatenated 'coeffData' to Improve Readability
% v2.0 - Updated to Support 'forceCoeffsExtended' Functionality


%% Case Selection

disp('Case Selection');
disp('---------------');

caseFolder = uigetdir('~/OpenFOAM/', 'Select Case');

campaignID = caseFolder((strfind(caseFolder, 'results/') + 8):(max(strfind(caseFolder, '/')) - 1));
caseID = caseFolder((max(strfind(caseFolder, '/')) + 1):end);

disp(' ');

disp(['Case: ', caseID]);

disp(' ');
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
    else
        
        % Sort Time Directories
        [~, index] = sort(str2double({coeffDirs.name}));
        coeffDirs = coeffDirs(index);
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


%% Perform Blockage Correction

if blockageCorrection
    
    if strcmp(campaignID, 'Windsor_fullScale')
        U = 22.222222222222222; % m/s
        rho = 1.269; % kg/m^3
        
        Am = ((4 * 0.289) * (4 * 0.389)) + (2 * (((4 * 0.05) - 0.018) * (4 * 0.055)));
        At = 14.336 * 26.624;
    elseif strcmp(campaignID, 'Windsor_Upstream_2023')
        U = 40; % m/s
        rho = 1.269; % kg/m^3
        
        Am = (0.289 * 0.389) + (2 * (0.046 * 0.055));
        At = (2 * (0.9519083 + (3.283 * tan(atan(0.0262223 / 9.44)))) * 1.32);
    end
    
    Fl = coeffData.Cl * 0.5 * rho * U^2 * Am;
    Fd = coeffData.Cd * 0.5 * rho * U^2 * Am;
    Fs = coeffData.Cs * 0.5 * rho * U^2 * Am;
    
    Ucorr = U * (At / (At - Am));
    
    coeffData.Cl = Fl / (0.5 * rho * Ucorr^2 * Am); clear Fl;
    coeffData.Cd = Fd / (0.5 * rho * Ucorr^2 * Am); clear Fd;
    coeffData.Cs = Fs / (0.5 * rho * Ucorr^2 * Am); clear Fs;
end


%% Calculate Time Average

Cl_Mean = zeros([height(coeffData.time),1]);
Cd_Mean = Cl_Mean;
Cs_Mean = Cl_Mean;

if strcmp(campaignID, 'Windsor_fullScale')
    startTime = 8;
elseif strcmp(campaignID, 'Windsor_Upstream_2023')
    startTime = 1;
end

startSample = find(coeffData.time == startTime);

for i = 1:height(coeffData.time)
    
    if i < startSample
        Cl_Mean(i) = nan;
        Cd_Mean(i) = nan;
        Cs_Mean(i) = nan;
    else
        Cl_Mean(i) = mean(coeffData.Cl(startSample:i));
        Cd_Mean(i) = mean(coeffData.Cd(startSample:i));
        Cs_Mean(i) = mean(coeffData.Cs(startSample:i));
    end
    
end


%% Plot Drag

% Initialise Figure
fig = fig + 1;
figName = 'Drag_Coefficient_Convergence';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 20, 'layer', 'top');

% Plot Coefficient
plot(coeffData.time, coeffData.Cd, 'color', graphColours(1), 'lineWidth', 2);
plot(coeffData.time, Cd_Mean, 'color', graphColours(2), 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;

if strcmp(campaignID, 'Windsor_fullScale')
    xlim([8; 32]);
%     ylim([0.26; 0.36]);
    tickData = (12.8:4.8:27.2);
    xticks(tickData);
%     tickData = (0.28:0.02:0.34);
%     yticks(tickData);
    xtickformat('%.1f');
%     ytickformat('%.2f');
elseif strcmp(campaignID, 'Windsor_Upstream_2023')
    xlim([1; 4]);
    ylim([0.3; 0.4]);
    tickData = (1.6:0.6:3.4);
    xticks(tickData);
    tickData = (0.32:0.02:0.38);
    yticks(tickData);
    xtickformat('%.1f');
    ytickformat('%.2f');
end

xlabel({'{Time $(s)$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$C_{_{D}}$}'}, 'interpreter', 'latex');
legend({'Instantaneous', ...
        'Expanding Average'}, ...
       'location', 'northWest', 'orientation', 'vertical', 'interpreter', 'latex', ...
       'fontSize', 18, 'box', 'off');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');


%% Plot Side Force

% Initialise Figure
fig = fig + 1;
figName = 'Side_Coefficient';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 20, 'layer', 'top');

% Plot Coefficient
plot(coeffData.time, coeffData.Cs, 'color', graphColours(1), 'lineWidth', 2);
% plot(coeffData.time, Cs_Mean, 'color', graphColours(2), 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;

if strcmp(campaignID, 'Windsor_fullScale')
    xlim([8; 32]);
    ylim([-0.06; 0.06]);
    tickData = (12.8:4.8:27.2);
    xticks(tickData);
    tickData = (-0.036:0.024:0.036);
    yticks(tickData);
    xtickformat('%.1f');
    ytickformat('%+.3f');
elseif strcmp(campaignID, 'Windsor_Upstream_2023')
    xlim([1; 4]);
    ylim([-0.06; 0.06]);
    tickData = (12.8:4.8:27.2);
    xticks(tickData);
    tickData = (-0.036:0.024:0.036);
    yticks(tickData);
    xtickformat('%.1f');
    ytickformat('%+.3f');
end

xlabel({'{Time $(s)$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$C_{_{S}}$}'}, 'interpreter', 'latex');
% legend({'\textit{Instantaneous}', ...
%         '\textit{Expanding Average}'}, ...
%        'location', 'northWest', 'orientation', 'vertical', 'interpreter', 'latex', ...
%        'fontSize', 16, 'box', 'off');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');