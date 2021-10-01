%% Lagrangian Mass Flux Planar POD Calculator v1.0

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('===============================================');
disp ('Lagrangian Mass Flux Planar POD Calculator v1.0');
disp ('===============================================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit


%% Case Initialisation

[caseFolder, xDims, yDims, zDims, timeDirs, deltaT, geometry] = initialiseCase('global');

disp(' ');
disp(' ');


%% Mapping Method

disp('MAPPING LOCATION');
disp('----------------');

disp(' ');
disp('Possible Mapping Locations:');
disp('    A: Base');
disp('    B: Far-Field Extraction Plane');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Mapping Location [A/B]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        method = 'A';
        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        method = 'B';
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end

disp(' ');
disp(' ');


%% Map Data Acquisition

disp('CONTAMINANT MAP ACQUISITION');
disp('---------------------------');

valid = false;
while ~valid
    disp(' ');
    [fileName, filePath] = uigetfile('/mnt/Processing/Data/Numerical/MATLAB/mapDataInstantaneous/*.*', ...
                                     'Select Map Data');
                            
    if contains(filePath, '/MATLAB/mapDataInstantaneous/')
        
        switch method
            
            case 'A'
                
                if contains(fileName, 'Base')
                    disp(['    Loading: ', fileName]);
                    load([filePath, fileName]);
                    valid = true;
                else
                    disp('    WARNING: Base Map Not Available in Specified File');
                    clear fileName filePath;
                end
                
            case 'B'
                
                if contains(fileName, 'Planar')
                    disp(['    Loading: ', fileName]);
                    load([filePath, fileName]);
                    valid = true;
                else
                    disp('    WARNING: Planar Map Not Available in Specified File');
                    clear fileName filePath;
                end
                
        end
        
    else
        disp('    WARNING: Invalid File Selection');
        clear fileName filePath;
    end
    
end

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

disp('    Performing Planar POD Using the Snapshot Method')

% Shift Data Origin
if contains(caseFolder, 'Upstream')
    xDims = xDims + 1.325;
end

% Set and Normalise Dimensions
xPre = max(width(extractAfter(num2str(xDims(1), 8), '.')), width(extractAfter(num2str(xDims(2), 8), '.')));
yPre = max(width(extractAfter(num2str(yDims(1), 8), '.')), width(extractAfter(num2str(yDims(2), 8), '.')));
zPre = max(width(extractAfter(num2str(zDims(1), 8), '.')), width(extractAfter(num2str(zDims(2), 8), '.')));

if contains(caseFolder, 'Test_Block') || contains(caseFolder, 'Windsor')
    xDims = round(xDims / 1.044, xPre);
    yDims = round(yDims / 1.044, yPre);
    zDims = round(zDims / 1.044, zPre);
    
    switch method
        
        case 'A'
            xLims = round([-0.61075; 0.53325] / 1.044, xPre);
            yLims = round([-0.2445; 0.2445] / 1.044, yPre);
            zLims = round([0; 0.389] / 1.044, zPre);
            
        case 'B'
            planePosition = mapData.global.positionGrid(1,1);
            
            % Standard Limits
            xLims = [round(-0.96075 / 1.044, xPre); planePosition];
            yLims = round([-0.5945; 0.5945] / 1.044, yPre);
            zLims = round([0; 0.789] / 1.044, zPre);
            
%             % Tight for Velocity Probe Comparison
%             xLims = [round(-0.96075 / 1.044, xPre); planePosition];
%             yLims = round([-0.256; 0.256] / 1.044, yPre);
%             zLims = round([0.0005; 0.4005] / 1.044, zPre);
            
    end
    
else
    % Add Support for Future Geometries
end

% Identify Model Boundaries
switch method
    
    case 'A'
        part = fieldnames(geometry);
        for i = 1:height(part)

            if round(max(geometry.(part{i,1}).vertices(:,1)), xPre) == xDims(2)
                break
            end

            if i == height(part)
                disp(' ');
                error('Mismatch Between ''xDims'' and Geometry Bounding Box')
            end

        end

        geoPoints = unique(geometry.(part{i,1}).vertices, 'stable', 'rows');
        index = find(round(geoPoints(:,1), xPre) == xDims(2));

        xDimsBase = xDims(2);
        yDimsBase = round([min(geoPoints(index,2)); max(geoPoints(index,2))], yPre);
        zDimsBase = round([min(geoPoints(index,3)); max(geoPoints(index,3))], zPre);
        
        yDimsBase(1) = yDimsBase(1) + 0.002;
        yDimsBase(2) = yDimsBase(2) - 0.002;
        zDimsBase(1) = zDimsBase(1) + 0.002;
        zDimsBase(2) = zDimsBase(2) - 0.002;
        
        modelOutline = [];
        
    case 'B'
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

        modelOutline(:,2) = round(modelOutline(:,2) / 1.044, yPre);
        modelOutline(:,3) = round(modelOutline(:,3) / 1.044, zPre);
        
end

% Initialise POD Variables
PODdata.time = mapData.global.time(:,1);
PODdata.position = mapData.global.positionGrid;

% Calculate Time-Averaged Deposition
PODdata.massMean = zeros(height(PODdata.position),1);

for i = 1:height(mapData.global.time)
    PODdata.massMean = PODdata.massMean + mapData.global.mass{i,1};
end

PODdata.massMean = PODdata.massMean / height(PODdata.time);

PODdata.massMeanNorm = PODdata.massMean / 2.933392467046334e-07; % Square-Back Base Cumulative Peak

% Calculate CoM of Time-Averaged Deposition
PODdata.massMeanCoM = zeros(1,3);

switch method
    
    case 'A'
        PODdata.massMeanCoM(1,1) = xDimsBase + 0.002;
        PODdata.massMeanCoM(1,2) = sum(PODdata.massMean .* PODdata.position(:,2)) / ...
                                                sum(PODdata.massMean);
        PODdata.massMeanCoM(1,3) = sum(PODdata.massMean .* PODdata.position(:,3)) / ...
                                                sum(PODdata.massMean);
        
    case 'B'
        PODdata.massMeanCoM(1,1) = planePosition;
        PODdata.massMeanCoM(1,2) = sum(PODdata.massMean .* PODdata.position(:,2)) / ...
                                       sum(PODdata.massMean);
        PODdata.massMeanCoM(1,3) = sum(PODdata.massMean .* PODdata.position(:,3)) / ...
                                       sum(PODdata.massMean);
                                            
end

% Calculate Calculate Instantaneous Deposition Fluctuations
for i = 1:height(mapData.global.time)
    PODdata.massPrime{i,1} = mapData.global.mass{i,1} - PODdata.massMean;
end

% Assemble Snapshot Matrix
Ns = height(PODdata.position); % Number of Spatial Points
Nt = height(PODdata.time); % Number of Temporal Instances

PODdata.snapshotMatrix = zeros(Nt,Ns);

for i = 1:Nt
    
    for j = 1:Ns
        PODdata.snapshotMatrix(i,j) = PODdata.massPrime{i,1}(j,1);
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

% Clear Redundant Data
clearvars mapData;

disp(' ');


%% Presenting Mode Energy Content

disp('    Presenting Mode Energy Content')
disp(['        First ', num2str(modesEnergetic), ' Modes Each Contain Greater Than 1% of Total Energy']);
disp(['        First ', num2str(modes80percent), ' Modes Contain Approximately 80% of Total Energy']);

% Figure Setup
fig = fig + 1;

switch method
    
    case 'A'
        figName = 'Mass_Flux_POD_Base_Mode_Energy_Content';
        set(figure(fig), 'outerPosition', [25, 25, 850, 850], 'name', figName);
        set(gca, 'fontName', 'LM Mono 12', 'fontSize', 18, 'layer', 'top');
        hold on;
        
    case 'B'
        figName = 'Mass_Flux_POD_Planar_Mode_Energy_Content';
        set(figure(fig), 'outerPosition', [25, 25, 850, 850], 'name', figName);
        set(gca, 'fontName', 'LM Mono 12', 'fontSize', 18, 'layer', 'top');
        hold on;
        
end

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
disp(' ');

tic;

if plotModes == 0
    disp('  Skipping Mode Presentation');
else
    % Identify Individual Modes    
    for i = plotModes
        disp(['    Presenting Mode #', num2str(i)])
        mode = ['M', num2str(i)];
        
        switch method
            
            case 'A'
                xLimsData = xDimsBase + 0.002;
                yLimsData = yDimsBase;
                zLimsData = zDimsBase;
                positionData = PODdata.position;
                depositionData = PODdata.phi_mode(:,i);
                figName = ['Mass_Flux_POD_Base_M', num2str(i)];
                cMap = turbo(24);
                CoM = [];
                
                figTitle = {' ', ' '};
%                 figTitle = {'Energy Content', [num2str(PODdata.modeEnergy(i,1)), '%']};
            
                cLims = 'auto';
                xLimsPlot = xLims;
                yLimsPlot = yLims;
                zLimsPlot = zLims;

                fig = contaminantPlots(xLimsData, yLimsData, zLimsData, ...
                                       positionData, depositionData, ...
                                       fig, figName, cMap, geometry, ...
                                       CoM, modelOutline, figTitle, cLims, ...
                                       xLimsPlot, yLimsPlot, zLimsPlot);
            
            case 'B'
                xLimsData = planePosition;
                yLimsData = yLims;
                zLimsData = zLims;
                positionData = PODdata.position;
                depositionData = PODdata.phi_mode(:,i);
                figName = ['Mass_Flux_POD_Planar_M', num2str(i)];
                cMap = turbo(24);
                CoM = [];
                
                figTitle = {' ', ' '};
%                 figTitle = {'Energy Content', [num2str(PODdata.modeEnergy(i,1)), '%']};

                cLims = 'auto';
                xLimsPlot = xLims;
                yLimsPlot = yLims;
                zLimsPlot = zLims;

                fig = contaminantPlots(xLimsData, yLimsData, zLimsData, ...
                                       positionData, depositionData, ...
                                       fig, figName, cMap, geometry, ...
                                       CoM, modelOutline, figTitle, cLims, ...
                                       xLimsPlot, yLimsPlot, zLimsPlot);
                
        end

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
   

%% Deposition Reconstruction

disp('DEPOSITION RECONSTRUCTION');
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
        plotModes = 1:Nt;
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
        disp(' ');
    end

end

disp(' ');

disp('***********');
disp('  Running  ');
disp(' ');

tic;

% Present Mean Deposition
disp('    Presenting Mean Deposition Map');

switch method
    
    case 'A'
        xLimsData = xDimsBase + 0.002;
        yLimsData = yDimsBase;
        zLimsData = zDimsBase;
        positionData = PODdata.position;
        depositionData = PODdata.massMeanNorm;
        figName = 'Mass_Flux_POD_Base_Mean';
        cMap = viridis(24);
        CoM = PODdata.massMeanCoM;
        figTitle = {' ', ' '};
        cLims = [0, 0.002]; % Windsor Full Duration
        xLimsPlot = xLims;
        yLimsPlot = yLims;
        zLimsPlot = zLims;

        fig = contaminantPlots(xLimsData, yLimsData, zLimsData, ...
                               positionData, depositionData, ...
                               fig, figName, cMap, geometry, ...
                               CoM, modelOutline, figTitle, cLims, ...
                               xLimsPlot, yLimsPlot, zLimsPlot);

    case 'B'
        xLimsData = planePosition;
        yLimsData = yLims;
        zLimsData = zLims;
        positionData = PODdata.position;
        depositionData = PODdata.massMeanNorm;
        figName = 'Mass_Flux_POD_Planar_Mean';
        cMap = viridis(24);
        CoM = PODdata.massMeanCoM;
        figTitle = {' ', ' '};
        cLims = [0, 0.02]; % Windsor Full Duration
        xLimsPlot = xLims;
        yLimsPlot = yLims;
        zLimsPlot = zLims;

        fig = contaminantPlots(xLimsData, yLimsData, zLimsData, ...
                               positionData, depositionData, ...
                               fig, figName, cMap, geometry, ...
                               CoM, modelOutline, figTitle, cLims, ...
                               xLimsPlot, yLimsPlot, zLimsPlot);

end

disp(' ');

if plotModes == 0
    disp('    Skipping Mode Reconstruction');
    reconstruction = [];
else    
    % Perform Reconstruction
    reconstruction.time = PODdata.time;
    reconstruction.position = PODdata.position;
    
    for i = 1:Nt
        reconstruction.mass{i,1} = PODdata.massMean;
    end
    
    for i = plotModes
        disp(['    Adding Mode #', num2str(i), ' to Reconstruction']);
        
        % Identify Mode
        mode = ['M', num2str(i)];
        reconstruction.(mode).modeMatrix = PODdata.A_coeff(:,i) * PODdata.phi_mode(:,i)';
        
        for j = 1:Nt
        
            for k = 1:Ns
                reconstruction.(mode).massPrime{j,1}(k,1) = reconstruction.(mode).modeMatrix(j,k);
            end
            
        end
        
        % Add Mode to Reconstruction
        for j = 1:Nt
            reconstruction.mass{j,1} = reconstruction.mass{j,1} + ...
                                         reconstruction.(mode).massPrime{j,1};
        end
        
    end
        
    % Normalise and Calculate Reconstructed CoM
    for j = 1:Nt
%         reconstruction.massNorm{j,1} = reconstruction.mass{j,1} / max(cellfun(@max, reconstruction.mass));

        reconstruction.massNorm{j,1} = reconstruction.mass{j,1} / 2.933392467046334e-07; % Square-Back Base Cumulative Peak

        reconstruction.CoM{j,1} = zeros(1,3);
        
        switch method

            case 'A'
                reconstruction.CoM{j,1}(1,1) = xDimsBase + 0.002;
                reconstruction.CoM{j,1}(1,2) = sum(reconstruction.mass{j,1} .* reconstruction.position(:,2)) / ...
                                               sum(reconstruction.mass{j,1});
                reconstruction.CoM{j,1}(1,3) = sum(reconstruction.mass{j,1} .* reconstruction.position(:,3)) / ...
                                               sum(reconstruction.mass{j,1});

            case 'B'
                reconstruction.CoM{j,1}(1,1) = planePosition;
                reconstruction.CoM{j,1}(1,2) = sum(reconstruction.mass{j,1} .* reconstruction.position(:,2)) / ...
                                               sum(reconstruction.mass{j,1});
                reconstruction.CoM{j,1}(1,3) = sum(reconstruction.mass{j,1} .* reconstruction.position(:,3)) / ...
                                               sum(reconstruction.mass{j,1});

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

        disp(['        T = ', num2str(str2double(reconstruction.time{i,1}), '%.4f'), ' s']);

        switch method

            case 'A'
                xLimsData = xDimsBase + 0.002;
                yLimsData = yDimsBase;
                zLimsData = zDimsBase;
                positionData = reconstruction.position;
                depositionData = reconstruction.massNorm{i,1};
                fig = figHold;
                figTime = num2str(str2double(reconstruction.time{i,1}), '%.4f');
                figName = ['Mass_Flux_POD_Base_Reconstruction_T', ...
                           erase(figTime, '.')];
                cMap = viridis(24);
                CoM = reconstruction.CoM{i,1};
                
                figTitle = {' ', ' '};
%                 figTitle = {'Time (\it{s})', figTime};
            
                cLims = [0, max(cellfun(@max, reconstruction.massNorm))];
                xLimsPlot = xLims;
                yLimsPlot = yLims;
                zLimsPlot = zLims;

                fig = contaminantPlots(xLimsData, yLimsData, zLimsData, ...
                                       positionData, depositionData, ...
                                       fig, figName, cMap, geometry, ...
                                       CoM, modelOutline, figTitle, cLims, ...
                                       xLimsPlot, yLimsPlot, zLimsPlot);

            case 'B'
                xLimsData = planePosition;
                yLimsData = yLims;
                zLimsData = zLims;
                positionData = reconstruction.position;
                depositionData = reconstruction.massNorm{i,1};
                fig = figHold;
                figTime = num2str(str2double(reconstruction.time{i,1}), '%.4f');
                figName = ['Mass_Flux_POD_Planar_Reconstruction_T', ...
                           erase(figTime, '.')];
                cMap = viridis(24);
                CoM = reconstruction.CoM{i,1};
                
                figTitle = {' ', ' '};
%                 figTitle = {'Time (\it{s})', figTime};
            
                cLims = [0, max(cellfun(@max, reconstruction.massNorm))];
                xLimsPlot = xLims;
                yLimsPlot = yLims;
                zLimsPlot = zLims;

                fig = contaminantPlots(xLimsData, yLimsData, zLimsData, ...
                                       positionData, depositionData, ...
                                       fig, figName, cMap, geometry, ...
                                       CoM, modelOutline, figTitle, cLims, ...
                                       xLimsPlot, yLimsPlot, zLimsPlot);

        end

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

        if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/massFluxPlanarPODdata/', caseFolder(namePos:end)], 'dir')
            mkdir(['/mnt/Processing/Data/Numerical/MATLAB/massFluxPlanarPODdata/', caseFolder(namePos:end)]);
        end

        switch method

            case 'A'
                mapType = 'Base';

            case 'B'
                mapType = 'Planar';

        end

        startInst = erase(num2str(str2double(PODdata.time{1,1}), '%.4f'), '.');
        endInst = erase(num2str(str2double(PODdata.time{end,1}), '%.4f'), '.');
        
        disp(['    Saving to: /mnt/Processing/Data/Numerical/MATLAB/massFluxPlanarPODdata/', caseFolder(namePos:end), '/T', startInst, '_T', endInst, '_D', num2str(minD), '_D', num2str(maxD), '_' mapType, '.mat']);
        save(['/mnt/Processing/Data/Numerical/MATLAB/massFluxPlanarPODdata/', caseFolder(namePos:end), '/T', startInst, '_T', endInst, '_D', num2str(minD), '_D', num2str(maxD), '_' mapType, '.mat'], 'PODdata', 'minD', 'maxD', '-v7.3', '-noCompression');

        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end


%% Cleaning

clearvars -except PODdata reconstruction minD maxD
disp(' ');


%% Local Functions

function M = inputModes

    valid = false;
    while ~valid
        M = str2num(input('    List Desired Modes (Row Vector Form): ', 's')); %#ok<ST2NM>
        
        if any(isnan(M)) || ~isrow(M)
            disp('        WARNING: Invalid Entry');
        else
            valid = true;
        end
        
    end
    
end