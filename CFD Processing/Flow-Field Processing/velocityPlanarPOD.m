%% Velocity Planar POD Calculator v1.1

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('===================================');
disp ('Velocity Planar POD Calculator v1.1');
disp ('===================================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Formatting Changes


%% Case Initialisation

[caseFolder, xDims, yDims, zDims, timeDirs, deltaT, geometry] = initialiseCase('PODprobe');

disp(' ');
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
        [probeData] = PODprobeData(caseFolder, timeDirs);
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        [fileName, filePath] = uigetfile('/mnt/Processing/Data/Numerical/MATLAB/PODprobeData/*.*', 'Select Probe Data');

        if contains(filePath, '/MATLAB/PODprobeData/')
            disp(['    Loading: ', fileName]);
            load([filePath, fileName])
            valid = true;
        else
            disp('    WARNING: Invalid File Selection');
            clearvars fileName filePath;
        end

    else
        disp('    WARNING: Invalid Entry');
    end

end

if contains(caseFolder, 'Windsor')

    % Shift Data Origin
    if contains(caseFolder, 'Upstream')
        probeData.position(:,1) = probeData.position(:,1) + 1.325;
    else
        % Add Support for Future Configurations
    end
    
else
    % Add Support for Future Geometries
end

% Identify Data Limits
xLimsData = [min(probeData.position(:,1)); max(probeData.position(:,1))];
yLimsData = [min(probeData.position(:,2)); max(probeData.position(:,2))];
zLimsData = [min(probeData.position(:,3)); max(probeData.position(:,3))];
           
disp(' ');


%% Specify Plane

disp('SPECIFY PLANE');
disp('-------------');

disp(' ');

disp('Probe Volume:');
disp(['    X: ', num2str(xLimsData(1)), ' [m] -> ', num2str(xLimsData(2)), ' [m]']);
disp(['    Y: ', num2str(yLimsData(1)), ' [m] -> ', num2str(yLimsData(2)), ' [m]']);
disp(['    Z: ', num2str(zLimsData(1)), ' [m] -> ', num2str(zLimsData(2)), ' [m]']);

disp(' ');

disp('Possible Plane Orientations:');
disp('    X: Normal [1 0 0]');
disp('    Y: Normal [0 1 0]');
disp('    Z: Normal [0 0 1]');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Plane Orientation [X/Y/Z]: ', 's');

    if selection == 'x' | selection == 'X' %#ok<OR2>
        planeOrientation = 'X';
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        planeOrientation = 'Y';
        valid = true;
    elseif selection == 'z' | selection == 'Z' %#ok<OR2>
        planeOrientation = 'Z';
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end

xLimsPlane = [0; 0];
yLimsPlane = [0; 0];
zLimsPlane = [0; 0];

valid = false;
while ~valid
    
    switch planeOrientation
        
        case 'X'
            xLimsPlane(1) = inputPos('Planar', 'x');
            xLimsPlane(2) = xLimsPlane(1);

            yLimsPlane(1) = inputPos('Lower', 'y');
            yLimsPlane(2) = inputPos('Upper', 'y');
            yLimsPlane = sort(yLimsPlane);

            zLimsPlane(1) = inputPos('Lower', 'z');
            zLimsPlane(2) = inputPos('Upper', 'z');
            zLimsPlane = sort(zLimsPlane);

        case 'Y'
            xLimsPlane(1) = inputPos('Lower', 'x');
            xLimsPlane(2) = inputPos('Upper', 'x');
            xLimsPlane = sort(xLimsPlane);

            yLimsPlane(1) = inputPos('Planar', 'y');
            yLimsPlane(2) = yLimsPlane(1);

            zLimsPlane(1) = inputPos('Lower', 'z');
            zLimsPlane(2) = inputPos('Upper', 'z');
            zLimsPlane = sort(zLimsPlane);

        case 'Z'
            xLimsPlane(1) = inputPos('Lower', 'x');
            xLimsPlane(2) = inputPos('Upper', 'x');
            xLimsPlane = sort(xLimsPlane);

            yLimsPlane(1) = inputPos('Lower', 'y');
            yLimsPlane(2) = inputPos('Upper', 'y');
            yLimsPlane = sort(yLimsPlane);

            zLimsPlane(1) = inputPos('Planar', 'z');
            zLimsPlane(2) = zLimsPlane(1);
            
    end
    
    if xLimsPlane(1) < xLimsData(1) || xLimsPlane(2) > xLimsData(2) || ...
       yLimsPlane(1) < yLimsData(1) || yLimsPlane(2) > yLimsData(2) || ...
       zLimsPlane(1) < zLimsData(1) || zLimsPlane(2) > zLimsData(2)
       disp('        WARNING: Plane Lies Outside Probe Volume');
       disp(' ');
    else
        valid = true;
    end
    
end


%% Identify Planar Probe Data

% Shift Requested Plane to Nearest Probe Plane
switch planeOrientation
    
    case 'X'
        planePositions = unique(probeData.position(:,1), 'stable');
        [offset, index] = min(abs(planePositions - xLimsPlane(1)));
        
        if offset ~= 0
            disp(' ');
            disp('WARNING: Requested Plane Unavailable');
            disp(['    Shifting X: ', num2str(xLimsPlane(1)), ' [m] -> ', num2str(planePositions(index)), ' [m]']);
            xLimsPlane(1) = planePositions(index);
            xLimsPlane(2) = xLimsPlane(1);
        end
        
        planePosition = xLimsPlane(1);
            
        index = find(probeData.position(:,1) == xLimsPlane(1) & ...
                     probeData.position(:,2) >= yLimsPlane(1) & probeData.position(:,2) <= yLimsPlane(2) & ...
                     probeData.position(:,3) >= zLimsPlane(1) & probeData.position(:,3) <= zLimsPlane(2));
        
    case 'Y'
        planePositions = unique(probeData.position(:,2), 'stable');
        [offset, index] = min(abs(planePositions - yLimsPlane(1)));
        
        if offset ~= 0
            disp(' ');
            disp('WARNING: Requested Plane Unavailable');
            disp(['    Shifting Y: ', num2str(yLimsPlane(1)), ' [m] -> ', num2str(planePositions(index)), ' [m]']);
            yLimsPlane(1) = planePositions(index);
            yLimsPlane(2) = yLimsPlane(1);
        end
            
        planePosition = yLimsPlane(1);
        
        index = find(probeData.position(:,1) >= xLimsPlane(1) & probeData.position(:,1) <= xLimsPlane(2) & ...
                     probeData.position(:,2) == yLimsPlane(1) & ...
                     probeData.position(:,3) >= zLimsPlane(1) & probeData.position(:,3) <= zLimsPlane(2));
        
    case 'Z'
        planePositions = unique(probeData.position(:,3), 'stable');
        [offset, index] = min(abs(planePositions - zLimsPlane(1)));
        
        if offset ~= 0
            disp(' ');
            disp('WARNING: Requested Plane Unavailable');
            disp(['    Shifting Z: ', num2str(zLimsPlane(1)), ' [m] -> ', num2str(planePositions(index)), ' [m]']);
            zLimsPlane(1) = planePositions(index);
            zLimsPlane(2) = zLimsPlane(1);
        end
        
        planePosition = zLimsPlane(1);
        
        index = find(probeData.position(:,1) >= xLimsPlane(1) & probeData.position(:,1) <= xLimsPlane(2) & ...
                     probeData.position(:,2) >= yLimsPlane(1) & probeData.position(:,2) <= yLimsPlane(2) & ...
                     probeData.position(:,3) == zLimsPlane(1));
        
end

% Identify Model Boundaries
if contains(caseFolder, 'Windsor')
    
    switch planeOrientation
        
        case 'X'
            
            if (planePosition < xDims(1) || planePosition > xDims(2)) && ...
               (yLimsPlane(1) < yDims(1) && yLimsPlane(2) > yDims(2)) && ...
               (zLimsPlane(1) < zDims(1) && zLimsPlane(2) > zDims(2))
           
                modelOutline = [
                                planePosition, -0.1945, 0.004;
                                planePosition, -0.1945, 0.339;
                                planePosition, 0.1945, 0.339;
                                planePosition, 0.1945, 0.004;
                                planePosition, 0.1395, 0.004;
                                planePosition, 0.1395, 0.05;
                                planePosition, -0.1395 0.05;
                                planePosition, -0.1395, 0.004;
                                planePosition, -0.1945, 0.004;
                               ];
            else
                modelOutline = [];
            end
            
        case {'Y', 'Z'}
            modelOutline = [];
    
    end
    
else
     % Add Support for Future Geometries
end       

% Collect Planar Velocity Data
PODdata.time = probeData.time;
PODdata.position = probeData.position(index,:);
PODdata.u = [];
PODdata.v = [];
PODdata.w = [];
PODdata.uMean = probeData.uMean(index,:);
PODdata.vMean = probeData.vMean(index,:);
PODdata.wMean = probeData.wMean(index,:);
PODdata.uPrime = [];
PODdata.vPrime = [];
PODdata.wPrime = [];


for i = 1:height(probeData.time)
    PODdata.u{i,1} = probeData.u{i,1}(index,:);
    PODdata.v{i,1} = probeData.v{i,1}(index,:);
    PODdata.w{i,1} = probeData.w{i,1}(index,:);
    PODdata.uPrime{i,1} = probeData.uPrime{i,1}(index,:);
    PODdata.vPrime{i,1} = probeData.vPrime{i,1}(index,:);
    PODdata.wPrime{i,1} = probeData.wPrime{i,1}(index,:);
end

% Set and Normalise Dimensions
xPre = max(width(extractAfter(num2str(xDims(1), 8), '.')), width(extractAfter(num2str(xDims(2), 8), '.')));
yPre = max(width(extractAfter(num2str(yDims(1), 8), '.')), width(extractAfter(num2str(yDims(2), 8), '.')));
zPre = max(width(extractAfter(num2str(zDims(1), 8), '.')), width(extractAfter(num2str(zDims(2), 8), '.')));

if contains(caseFolder, 'Windsor')
    xLimsPlot = round([0.31875; max(xLimsPlane)] / 1.044, xPre);
    yLimsPlot = round([-0.2445; 0.2445] / 1.044, zPre);
    zLimsPlot = round([0; 0.389] / 1.044, zPre);
    
    xLimsPlot = round([0.31875; max(xLimsPlane)] / 1.044, xPre);
    yLimsPlot = round([-0.256; 0.256] / 1.044, zPre);
    zLimsPlot = round([0.0005; 0.4005] / 1.044, zPre);
else
    % Add Support for Future Geometries
end

PODdata.position(:,1) = round(PODdata.position(:,1) / 1.044, xPre);
PODdata.position(:,2) = round(PODdata.position(:,2) / 1.044, xPre);
PODdata.position(:,3) = round(PODdata.position(:,3) / 1.044, xPre);

xDims = round(xDims / 1.044, xPre);
yDims = round(yDims / 1.044, yPre);
zDims = round(zDims / 1.044, zPre);

if ~isempty(modelOutline)
    modelOutline(:,1) = round(modelOutline(:,1) / 1.044, xPre);
    modelOutline(:,2) = round(modelOutline(:,2) / 1.044, yPre);
    modelOutline(:,3) = round(modelOutline(:,3) / 1.044, zPre);
end

switch planeOrientation
    
    case 'X'
        planePosition = round(planePosition / 1.044, xPre);
        xLimsPlane = round(xLimsPlane / 1.044, xPre);
        yLimsPlane = round(yLimsPlane / 1.044, yPre);
        zLimsPlane = round(zLimsPlane / 1.044, zPre);
        
    case 'Y'
        planePosition = round(planePosition / 1.044, yPre);
        xLimsPlane = round(xLimsPlane / 1.044, xPre);
        yLimsPlane = round(yLimsPlane / 1.044, yPre);
        zLimsPlane = round(zLimsPlane / 1.044, zPre);
        
    case 'Z'
        planePosition = round(planePosition / 1.044, zPre);
        xLimsPlane = round(xLimsPlane / 1.044, xPre);
        yLimsPlane = round(yLimsPlane / 1.044, yPre);
        zLimsPlane = round(zLimsPlane / 1.044, zPre);
        
end

% Clear Redundant Data
clearvars probeData;

disp(' ');
disp(' ');


%% Planar POD Calculation (Snapshot Method)

disp('PLANAR POD');
disp('----------');

disp(' ');

disp('***********');
disp('  Running  ');

tic;

disp(' ');

% Assemble Snapshot Matrix
disp('    Performing Planar POD Using the Snapshot Method')

Ns = height(PODdata.position); % Number of Spatial Points
Nt = height(PODdata.time); % Number of Temporal Instances

PODdata.snapshotMatrix = zeros(Nt,(3 * Ns)); % U, V, W

for i = 1:Nt
    
    for j = 1:Ns
        PODdata.snapshotMatrix(i,j) = PODdata.uPrime{i,1}(j,1);
        PODdata.snapshotMatrix(i,(j + Ns)) = PODdata.vPrime{i,1}(j,1);
        PODdata.snapshotMatrix(i,(j + (2 * Ns))) = PODdata.wPrime{i,1}(j,1);
    end
    
end

% Produce Correlation Matrix
PODdata.C = (PODdata.snapshotMatrix * PODdata.snapshotMatrix') / (Nt - 1);

% Solve Eigenvalue Problem
[PODdata.A_mode, PODdata.lambda] = eig(PODdata.C, 'vector');

% Sort Eigenvalues and Eigenvectors in Descending Order
[PODdata.lambda, index] = sort(PODdata.lambda, 'descend');
PODdata.A_mode = PODdata.A_mode(:,index); % Temporal Modes

% Calculate Spatial Coefficients
PODdata.phi_coeff = PODdata.snapshotMatrix' * PODdata.A_mode;

% Normalisation to Match Direct Method
PODdata.phi_mode = normc(PODdata.phi_coeff); % Spatial Modes
PODdata.A_coeff = PODdata.snapshotMatrix * PODdata.phi_mode; % Temporal Coefficients

% Identify Relative Mode Energy Content
PODdata.modeEnergy = (PODdata.lambda / sum(PODdata.lambda)) * 100;
modesEnergetic = height(PODdata.modeEnergy(PODdata.modeEnergy > 1));
modes80percent = find(cumsum(PODdata.modeEnergy) > 80, 1);

disp(' ');


%% Mode Energy Content

disp('    Presenting Mode Energy Content')
disp(['        First ', num2str(modesEnergetic), ' Modes Each Contain Greater Than 1% of Total Energy']);
disp(['        First ', num2str(modes80percent), ' Modes Contain Approximately 80% of Total Energy']);

% Figure Setup
fig = fig + 1;
figName = ['Velocity_POD_', planeOrientation, '_', erase(num2str(planePosition, '%.4f'), '.'), '_Mode_Energy_Content'];
set(figure(fig), 'outerPosition', [25, 25, 850, 850], 'name', figName);
set(gca, 'fontName', 'LM Mono 12', 'fontSize', 18, 'layer', 'top');
hold on;

% Plot
bar(PODdata.modeEnergy(1:20,1), ...
    'faceColor', [0.21176, 0.06667, 0.38824], 'lineStyle', 'none');

% Figure Formatting
title(' ');
box off;
xlim([0, 21]);
ylim([0, (ceil(PODdata.modeEnergy(1,1)) + 1)]);
tickData = 2:2:20;
xticks(tickData);
tickData = 0:(ceil(PODdata.modeEnergy(1,1)) + 1);
yticks(tickData(2:(end-1)));
xT = xlabel('Mode');
yT = ylabel('Energy Content (\it{%})');
xT.FontName = 'LM Roman 12';
yT.FontName = 'LM Roman 12';
hold off;

pause(2);
exportgraphics(gca, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);
% savefig(fig, ['~/MATLAB/Output/Figures/', figName]);

disp(' ');

executionTime = toc;

disp(['    Run Time: ', num2str(executionTime), 's']);
disp(' ');
disp('  Success  ')
disp('***********');

disp(' ');
disp(' ');


%% Mode Presentation

disp('MODE PRESENTATION');
disp('-----------------');

valid = false;
while ~valid
    disp(' ');
    selection = input(['Plot All ', num2str(modesEnergetic), ' Energetic Modes? [y/n]: '], 's');
    
    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotModes = inputModes;
        
        if (width(plotModes) > 1 && min(plotModes) <= 0) || max(plotModes) > Nt
            disp('        WARNING: Invalid Mode Selection');
        else
            valid = true;
        end
        
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotModes = 1:modesEnergetic;
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
        disp(' ');
    end

end

disp(' ');

disp('***********');
disp('  Running  ');

tic;

disp(' ');

if plotModes == 0
    disp('  Skipping Mode Presentation');
else    
    % Present Individual Modes
    for i = plotModes
        disp(['    Presenting Mode #', num2str(i)])
        mode = ['M', num2str(i)];
        
        % Single-Component Plots
        xLimsData = xLimsPlane;
        yLimsData = yLimsPlane;
        zLimsData = zLimsPlane;
        positionData = PODdata.position;
        velocityData = [PODdata.phi_mode((1:Ns),i), ...
                        PODdata.phi_mode((Ns + 1):(2 * Ns),i), ...
                        PODdata.phi_mode(((2* Ns) + 1):(3 * Ns),i)];
        nComponents = 1;
        
        for j = ['u', 'v', 'w']
            component = j;
            figName = ['Velocity_POD_', planeOrientation, '_', erase(num2str(planePosition, '%.4f'), '.'), ...
                       '_', mode, '_', component, 'Prime'];
            cMap = turbo(24);
            streamlines = false;

            figTitle = {' ', ' '};
    %         figTitle = {'Energy Content', [num2str(PODdata.modeEnergy(i,1)), '%']};

            cLims = 'auto';

            fig = velocityPlots(planeOrientation, planePosition, ...
                                 xLimsData, yLimsData, zLimsData, ...
                                 positionData, velocityData, ...
                                 nComponents, component, ...
                                 fig, figName, cMap, geometry, ...
                                 streamlines, modelOutline, figTitle, cLims, ...
                                 xLimsPlot, yLimsPlot, zLimsPlot);
        end
        
        % Three-Component Plots
        xLimsData = xLimsPlane;
        yLimsData = yLimsPlane;
        zLimsData = zLimsPlane;
        positionData = PODdata.position;
        velocityData = [PODdata.phi_mode((1:Ns),i), ...
                        PODdata.phi_mode((Ns + 1):(2 * Ns),i), ...
                        PODdata.phi_mode(((2* Ns) + 1):(3 * Ns),i)];
        nComponents = 1;
        
        switch planeOrientation

            case 'X'
                component = 'u';

            case 'Y'
                component = 'v';

            case 'Z'
                component = 'w';

        end
        
        figName = ['Velocity_POD_', planeOrientation, '_', erase(num2str(planePosition, '%.4f'), '.'), ...
                   '_', mode];
        cMap = turbo(24);
        streamlines = true;
        
        figTitle = {' ', ' '};
%         figTitle = {'Energy Content', [num2str(PODdata.modeEnergy(i,1)), '%']};

        cLims = 'auto';

        fig = velocityPlots(planeOrientation, planePosition, ...
                             xLimsData, yLimsData, zLimsData, ...
                             positionData, velocityData, ...
                             nComponents, component, ...
                             fig, figName, cMap, geometry, ...
                             streamlines, modelOutline, figTitle, cLims, ...
                             xLimsPlot, yLimsPlot, zLimsPlot);
    end

end

disp(' ');

executionTime = toc;

disp(['    Run Time: ', num2str(executionTime), 's']);
disp(' ');
disp('  Success  ')
disp('***********');

disp(' ');
disp(' ');


%% Flow-Field Reconstruction

disp('FLOW-FIELD RECONSTRUCTION');
disp('-------------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input(['Reconstruct 80% of Modal Energy Content (', num2str(modes80percent), ' Modes)? [y/n]: '], 's');
    
    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotModes = inputModes;
        
        if (width(plotModes) > 1 && min(plotModes) <= 0) || max(plotModes) > Nt
            disp('        WARNING: Invalid Mode Selection');
        else
            valid = true;
        end
        
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotModes = 1:modes80percent;
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
        disp(' ');
    end

end

disp(' ');

disp('***********');
disp('  Running  ');

tic;

disp(' ');

% Present Mean Flow-Field
disp('    Presenting Mean Flow-Field');

xLimsData = xLimsPlane;
yLimsData = yLimsPlane;
zLimsData = zLimsPlane;
positionData = PODdata.position;
velocityData = [PODdata.uMean / 40, ...
                PODdata.vMean / 40, ...
                PODdata.wMean / 40];
nComponents = 3;
component = [];
figName = ['Velocity_POD_', planeOrientation, '_', erase(num2str(planePosition, '%.4f'), '.'), ...
           '_Mean_Velocity_Magnitude'];
cMap = viridis(24);
streamlines = true;
figTitle = {' ', ' '};
cLims = [0, 1];

fig = velocityPlots(planeOrientation, planePosition, ...
                    xLimsData, yLimsData, zLimsData, ...
                    positionData, velocityData, ...
                    nComponents, component, ...
                    fig, figName, cMap, geometry, ...
                    streamlines, modelOutline, figTitle, cLims, ...
                    xLimsPlot, yLimsPlot, zLimsPlot);
                 
disp(' ');

if plotModes == 0
    disp('    Skipping Mode Reconstruction');
    reconstruction = [];
else    
    % Perform Reconstruction
    reconstruction.time = PODdata.time;
    reconstruction.position = PODdata.position;
    
    for i = 1:Nt
        reconstruction.u{i,1} = PODdata.uMean;
        reconstruction.v{i,1} = PODdata.vMean;
        reconstruction.w{i,1} = PODdata.wMean;
    end
    
    for i = plotModes
        disp(['    Adding Mode #', num2str(i), ' to Reconstruction']);
        
        % Identify Mode
        mode = ['M', num2str(i)];
        reconstruction.(mode).modeMatrix = PODdata.A_coeff(:,i) * PODdata.phi_mode(:,i)';
        
        for j = 1:Nt
            
            for k = 1:Ns
                reconstruction.(mode).uPrime{j,1}(k,1) = reconstruction.(mode).modeMatrix(j,k);
                reconstruction.(mode).vPrime{j,1}(k,1) = reconstruction.(mode).modeMatrix(j,(k + Ns));
                reconstruction.(mode).wPrime{j,1}(k,1) = reconstruction.(mode).modeMatrix(j,(k + (2 * Ns)));
            end
            
        end
        
        % Add Mode to Reconstruction
        for j = 1:Nt
            reconstruction.u{j,1} = reconstruction.u{j,1} + reconstruction.(mode).uPrime{j,1};
            reconstruction.v{j,1} = reconstruction.v{j,1} + reconstruction.(mode).vPrime{j,1};
            reconstruction.w{j,1} = reconstruction.w{j,1} + reconstruction.(mode).wPrime{j,1};
        end
        
    end
    
    disp(' ');
    
    % Present Reconstruction
    disp('    Presenting Reconstruction');
    
    figHold = fig;
    
    for i = 1:Nt

        if i ~= 1
            clf(fig);
        end

        disp(['        T = ', num2str(reconstruction.time(i,1), '%.4f'), ' s']);
        
        xLimsData = xLimsPlane;
        yLimsData = yLimsPlane;
        zLimsData = zLimsPlane;
        positionData = reconstruction.position;
        velocityData = [reconstruction.u{i,1} / 40, ...
                        reconstruction.v{i,1} / 40, ...
                        reconstruction.w{i,1} / 40];
        nComponents = 3;
        component = [];
        fig = figHold;
        figTime = num2str(reconstruction.time(i,1), '%.4f');
        figName = ['Velocity_POD_', planeOrientation, '_', erase(num2str(planePosition, '%.4f'), '.'), ...
                   '_Reconstruction_T', erase(figTime, '.')];
        cMap = viridis(24);
        streamlines = true;
        
        figTitle = {' ', ' '};
%         figTitle = {'Time (\it{s})', figTime};
        
        cLims = [0, 1];

        fig = velocityPlots(planeOrientation, planePosition, ...
                            xLimsData, yLimsData, zLimsData, ...
                            positionData, velocityData, ...
                            nComponents, component, ...
                            fig, figName, cMap, geometry, ...
                            streamlines, modelOutline, figTitle, cLims, ...
                            xLimsPlot, yLimsPlot, zLimsPlot);
    end
    
end

disp(' ');

executionTime = toc;

disp(['    Run Time: ', num2str(executionTime), 's']);
disp(' ');
disp('  Success  ')
disp('***********');

disp(' ');


%% Saving POD Data

valid = false;
while ~valid
    disp(' ');
    selection = input('Save POD Data for Future Use? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        namePos = max(strfind(caseFolder, '/')) + 1;

        if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/velocityPlanarPODdata/', caseFolder(namePos:end)], 'dir')
            mkdir(['/mnt/Processing/Data/Numerical/MATLAB/velocityPlanarPODdata/', caseFolder(namePos:end)]);
        end
        
        startInst = erase(num2str(str2double(PODdata.time{1,1}), '%.4f'), '.');
        endInst = erase(num2str(str2double(PODdata.time{end,1}), '%.4f'), '.');
        
        disp(['    Saving to: /mnt/Processing/Data/Numerical/MATLAB/velocityPlanarPODdata/', caseFolder(namePos:end), '/T', startInst, '_T', endInst, '.mat']);
        save(['/mnt/Processing/Data/Numerical/MATLAB/velocityPlanarPODdata/', caseFolder(namePos:end), '/T', startInst, '_T', endInst, '.mat'], 'PODdata', '-v7.3', '-noCompression');

        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end


%% Cleaning

clearvars -except PODdata reconstruction
disp(' ');


%% Local Functions

function pos = inputPos(type, plane)

    valid = false;
    while ~valid
        pos = str2double(input(['    ', type, ' ', plane, '-position [m]: '], 's'));
        
        if isnan(pos) || length(pos) > 1
            disp('        WARNING: Invalid Entry');
        else
            valid = true;
        end

    end

end


function M = inputModes

    valid = false;
    while ~valid
        disp(' ');
        M = str2num(input('    List Desired Modes (Row Vector Form): ', 's')); %#ok<ST2NM>
        
        if any(isnan(M)) || ~isrow(M)
            disp('        WARNING: Invalid Entry');
        else
            valid = true;
        end
        
    end
    
end