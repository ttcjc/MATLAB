%% Planar POD Calculator v1.0

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('==========================');
disp ('Planar POD Calculator v1.0');
disp ('==========================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit


%% Case Initialisation

% [caseFolder, ~, ~, ~, timeDirs, deltaT, geometry] = initialiseCase('PODprobe');

caseFolder = '/home/lunet/ttcjc/OpenFOAM/ttcjc-7/results/Windsor_Square_wW_Balance_SC';
% format = 'PODprobe';
% [timeDirs, deltaT] = timeDirectories(caseFolder, format);
load('geometryTest.mat');


% disp(' ');
disp(' ');


%% Probe Data Acquisition

disp('PROBE DATA ACQUISITION');
disp('----------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Load Saved Probe Data? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        disp(' ');
        disp(' ');
        [probeData] = PODprobeData(caseFolder, timeDirs);
        valid = true;
        disp(' ');
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        [probeData, path] = uigetfile('~/Documents/Engineering/PhD/Data/Numerical/MATLAB/PODprobeData/*.*', 'Select Probe Data');

        if contains(path, '/MATLAB/PODprobeData/')
            disp(' ');
            disp(['    Loading: ', probeData]);
            load([path, probeData])
            valid = true;
            disp(' ');
        else
            clearvars probeData path;
            disp('    WARNING: Invalid File Selection');
            disp(' ');
        end

    else
        disp('    WARNING: Invalid Entry');
        disp(' ');
    end

end

xDims = [min(probeData.position(:,1)); max(probeData.position(:,1))];
yDims = [min(probeData.position(:,2)); max(probeData.position(:,2))];
zDims = [min(probeData.position(:,3)); max(probeData.position(:,3))];


disp(' ');


%% Specify Plane

disp('SPECIFY PLANE');
disp('-------------');

disp(' ');
disp('Probe Volume:');
disp(['    X: ', num2str(xDims(1)), ' [m] -> ', num2str(xDims(2)), ' [m]']);
disp(['    Y: ', num2str(yDims(1)), ' [m] -> ', num2str(yDims(2)), ' [m]']);
disp(['    Z: ', num2str(zDims(1)), ' [m] -> ', num2str(zDims(2)), ' [m]']);

disp(' ');
disp('Possible Plane Orientations:');
disp('    X: Normal [1 0 0]');
disp('    Y: Normal [0 1 0]');
disp('    Z: Normal [0 0 1]');
disp(' ');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Plane Orientation [X/Y/Z]: ', 's');

    if selection == 'x' | selection == 'X' %#ok<OR2>
        orientation = 'X';
        valid = true;
        disp(' ');
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        orientation = 'Y';
        valid = true;
        disp(' ');
    elseif selection == 'z' | selection == 'Z' %#ok<OR2>
        orientation = 'Z';
        valid = true;
        disp(' ');
    else
        disp('    WARNING: Invalid Entry');
        disp(' ');
    end

end

xDimsPlane = [0; 0];
yDimsPlane = [0; 0];
zDimsPlane = [0; 0];

valid = false;
while ~valid
    
    switch orientation
        case 'X'
            xDimsPlane(1) = inputPos('Planar', 'x');
            xDimsPlane(2) = xDimsPlane(1);
            xDimsPlane = sort(xDimsPlane);

            yDimsPlane(1) = inputPos('Lower', 'y');
            yDimsPlane(2) = inputPos('Upper', 'y');
            yDimsPlane = sort(yDimsPlane);

            zDimsPlane(1) = inputPos('Lower', 'z');
            zDimsPlane(2) = inputPos('Upper', 'z');
            zDimsPlane = sort(zDimsPlane);

        case 'Y'
            xDimsPlane(1) = inputPos('Lower', 'x');
            xDimsPlane(2) = inputPos('Upper', 'x');
            xDimsPlane = sort(xDimsPlane);

            yDimsPlane(1) = inputPos('Planar', 'y');
            yDimsPlane(2) = yDimsPlane(1);
            yDimsPlane = sort(yDimsPlane);

            zDimsPlane(1) = inputPos('Lower', 'z');
            zDimsPlane(2) = inputPos('Upper', 'z');
            zDimsPlane = sort(zDimsPlane);

        case 'Z'
            xDimsPlane(1) = inputPos('Lower', 'x');
            xDimsPlane(2) = inputPos('Upper', 'x');
            xDimsPlane = sort(xDimsPlane);

            yDimsPlane(1) = inputPos('Lower', 'y');
            yDimsPlane(2) = inputPos('Upper', 'y');
            yDimsPlane = sort(yDimsPlane);

            zDimsPlane(1) = inputPos('Planar', 'z');
            zDimsPlane(2) = zDimsPlane(1);
            zDimsPlane = sort(zDimsPlane);
    end
    
    if xDimsPlane(1) < xDims(1) || xDimsPlane(2) > xDims(2) || ...
       yDimsPlane(1) < yDims(1) || yDimsPlane(2) > yDims(2) || ...
       zDimsPlane(1) < zDims(1) || zDimsPlane(2) > zDims(2)
       disp('        WARNING: Plane Lies Outside Probe Volume');
       disp(' ');
    else
        valid = true;
        disp(' ');
    end
    
end


%% Identify Planar Probe Data

% Shift Requested Plane to Nearest Probe Plane
switch orientation
    case 'X'
        planePositions = unique(probeData.position(:,1), 'stable');
        [offset, index] = min(abs(planePositions - xDimsPlane(1)));
        
        if offset ~= 0
            disp('    WARNING: Requested Plane Unavailable');
            disp(['             Shifting X: ', num2str(xDimsPlane(1)), ' [m] -> ', num2str(planePositions(index)), ' [m]']);
            xDimsPlane(1) = planePositions(index);
            xDimsPlane(2) = xDimsPlane(1);
        end
            
            index = find(probeData.position(:,1) == xDimsPlane(1) & ...
                         probeData.position(:,2) >= yDimsPlane(1) & probeData.position(:,2) <= yDimsPlane(2) & ...
                         probeData.position(:,3) >= zDimsPlane(1) & probeData.position(:,3) <= zDimsPlane(2));
        
    case 'Y'
        planePositions = unique(probeData.position(:,2), 'stable');
        [offset, index] = min(abs(planePositions - yDimsPlane(1)));
        
        if offset ~= 0
            disp('    WARNING: Requested Plane Unavailable');
            disp(['             Shifting Y: ', num2str(yDimsPlane(1)), ' [m] -> ', num2str(planePositions(index)), ' [m]']);
            yDimsPlane(1) = planePositions(index);
            yDimsPlane(2) = yDimsPlane(1);
        end
            
            index = find(probeData.position(:,1) >= xDimsPlane(1) & probeData.position(:,1) <= xDimsPlane(2) & ...
                         probeData.position(:,2) == yDimsPlane(1) & ...
                         probeData.position(:,3) >= zDimsPlane(1) & probeData.position(:,3) <= zDimsPlane(2));
        
    case 'Z'
        planePositions = unique(probeData.position(:,3), 'stable');
        [offset, index] = min(abs(planePositions - zDimsPlane(1)));
        
        if offset ~= 0
            disp('    WARNING: Requested Plane Unavailable');
            disp(['             Shifting Z: ', num2str(zDimsPlane(1)), ' [m] -> ', num2str(planePositions(index)), ' [m]']);
            zDimsPlane(1) = planePositions(index);
            zDimsPlane(2) = zDimsPlane(1);
        end
        
            index = find(probeData.position(:,1) >= xDimsPlane(1) & probeData.position(:,1) <= xDimsPlane(2) & ...
                         probeData.position(:,2) >= yDimsPlane(1) & probeData.position(:,2) <= yDimsPlane(2) & ...
                         probeData.position(:,3) == zDimsPlane(1));
        
end

% Collect Planar Velocity Data
planarPODdata.time = probeData.time;
planarPODdata.position = probeData.position(index,:);
planarPODdata.u = [];
planarPODdata.v = [];
planarPODdata.w = [];
planarPODdata.uMean = probeData.uMean(index,:);
planarPODdata.vMean = probeData.vMean(index,:);
planarPODdata.wMean = probeData.wMean(index,:);
planarPODdata.uPrime = [];
planarPODdata.vPrime = [];
planarPODdata.wPrime = [];


for i = 1:size(probeData.time,1)
    planarPODdata.u{i,1} = probeData.u{i,1}(index,:);
    planarPODdata.v{i,1} = probeData.u{i,1}(index,:);
    planarPODdata.w{i,1} = probeData.u{i,1}(index,:);
    planarPODdata.uPrime{i,1} = probeData.uPrime{i,1}(index,:);
    planarPODdata.vPrime{i,1} = probeData.vPrime{i,1}(index,:);
    planarPODdata.wPrime{i,1} = probeData.wPrime{i,1}(index,:);
end

disp(' ');
disp(' ');


%% Planar POD Calculation (Snapshot Method)

tic;
disp('***********');
disp('  Running  ');
disp(' ');
disp('    Performing Planar POD Using the Snapshot Method')

% Assemble Fluctuation Matrix
planarPODdata.primeMatrix = zeros((3 * size(planarPODdata.position,1)), size(planarPODdata.time,1));

for i = 1:size(planarPODdata.position,1)
    
    for j = 1:size(planarPODdata.time,1)
        k = 3 * i;
        planarPODdata.primeMatrix(k-2,j) = planarPODdata.uPrime{j,1}(i);
        planarPODdata.primeMatrix((k-1),j) = planarPODdata.vPrime{j,1}(i);
        planarPODdata.primeMatrix(k,j) = planarPODdata.wPrime{j,1}(i);
    end
    
end

planarPODdata.R = planarPODdata.primeMatrix' * planarPODdata.primeMatrix; % Autocovariance matrix
[planarPODdata.eV, planarPODdata.D] = eig(planarPODdata.R); % (eV) eigenvectors, (D) eigenvalues in diagonal matrix
[planarPODdata.L, planarPODdata.I] = sort(diag(planarPODdata.D)); % sort eignevalues in ascending order

for i = 1:length(planarPODdata.D)
    planarPODdata.eValue(length(planarPODdata.D) + 1 - i) = planarPODdata.L(i); % Eigenvalues sorted in descending order
    planarPODdata.eVec(:,length(planarPODdata.D) + 1 - i) = planarPODdata.eV(:, planarPODdata.I(i)); % Eigenvectors sorted in the same order
end

planarPODdata.eValue(length(planarPODdata.eValue)) = 0; % last eigenvalue should be zero
planarPODdata.modeEnergy = planarPODdata.eValue / sum(planarPODdata.eValue); % relative energy associated with mode

i = 1;
while planarPODdata.modeEnergy(i) * 100 >= 1
    tmp = planarPODdata.primeMatrix * planarPODdata.eVec(:,i);
    planarPODdata.phi(:,i) = tmp / norm(tmp);
    i = i + 1;
end

disp(' ');
disp(['    First ', num2str(i-1), ' Modes Contain Greater Than 1% Of Total Energy']);

planarPODdata.An = (planarPODdata.phi' * planarPODdata.primeMatrix)';
executionTime = toc;

disp(' ');
disp(['    Execution Time: ', num2str(executionTime), 's']);
disp(' ');
disp('  Success  ');
disp('***********');
disp(' ');

valid = false;
while ~valid
    disp(' ');
    selection = input('Save POD Data for Future Use? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        namePos = max(strfind(caseFolder, '/')) + 1;
        
        if ~exist(['~/Documents/Engineering/PhD/Data/Numerical/MATLAB/planarPODdata/', caseFolder(namePos(end):end)], 'dir')
            mkdir(['~/Documents/Engineering/PhD/Data/Numerical/MATLAB/planarPODdata/', caseFolder(namePos(end):end)]);
        end
        
        startInst = erase(num2str(planarPODdata.time(1,1), '%.4f'), '.');
        endInst = erase(num2str(planarPODdata.time(end,1), '%.4f'), '.');
        x1 = erase(num2str(xDimsPlane(1), '%.4f'), '.');
        x2 = erase(num2str(xDimsPlane(2), '%.4f'), '.');
        y1 = erase(num2str(yDimsPlane(1), '%.4f'), '.');
        y2 = erase(num2str(yDimsPlane(2), '%.4f'), '.');
        z1 = erase(num2str(zDimsPlane(1), '%.4f'), '.');
        z2 = erase(num2str(zDimsPlane(2), '%.4f'), '.');
        disp(' ');
        disp(['    Saving to: ~/Documents/Engineering/PhD/Data/Numerical/MATLAB/planarPODdata/', caseFolder(namePos(end):end), '/S', startInst, '_E', endInst, '_x', x1, '_y', y1, '_z', z1, '_x', x2, '_y', y2, '_z', z2, '.mat']);
        save(['~/Documents/Engineering/PhD/Data/Numerical/MATLAB/planarPODdata/', caseFolder(namePos(end):end), '/S', startInst, '_E', endInst, '_x', x1, '_y', y1, '_z', z1, '_x', x2, '_y', y2, '_z', z2, '.mat'], 'planarPODdata', '-v7.3', '-noCompression');
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
        disp(' ');
    end

end

% Figure Setup
fig = fig + 1;
figure(fig);
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800], 'name', 'POD: Mode Energy Content');

% Plot
bar(planarPODdata.modeEnergy(1:i-1) * 100);

% Figure Formatting
xlabel({' ', 'Mode'});
ylabel({'Energy Content (\it{%})', ' '});
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
         'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;


%% Reconstruction



%% Local Functions

function pos = inputPos(type, plane)

    valid = false;
    while ~valid
        disp(' ');
        pos = str2double(input(['    ', type, ' ', plane, '-position [m]: '], 's'));
        
        if isnan(pos) || length(pos) > 1
            disp('        WARNING: Invalid Entry');
        else
            valid = true;
        end

    end

end
