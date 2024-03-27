%% Planar Spray POD Field Reconstructor v1.2
% ----
% Perform Low-Order Reconstructions Using the Output of 'planarSprayPOD'


%% Preamble

run preamble;

%#ok<*UNRCH>

normDensity = true; % Normalise Spray Density in Plots

normDims = true; % Normalise Spatial Dimensions in Plots

figSave = false; % Save .fig File(s)

disp('=========================================');
disp('Planar Spray POD Field Reconstructor v1.2');
disp('=========================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit (Functionality Separated From 'planarSprayPOD' To Improve Efficiency)
% v1.1 - Minor Update to Shift Preamble Into Separate Script
% v1.2 - Update To Correct Inconsistent Normalisation Throughout Repository


%% Select Mapping Location

disp('Mapping Location');
disp('-----------------');

disp(' ');

disp('Possible Mapping Locations:');
disp('    A: Base Contamination');
disp('    B: Far-Field Spray Transport');
disp('    C: Experimental Spray Transport');

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Select Mapping Location [A/B/C]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        format = 'A';
        
        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        format = 'B';
        
        valid = true;
    elseif selection == 'c' | selection == 'C' %#ok<OR2>
        format = 'C';
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Acquire Contaminant Map

disp('Spray POD Acquisition');
disp('----------------------');

switch format
    
    case {'A', 'B'}
        
        valid = false;
        while ~valid
            disp(' ');
            
            [fileName, filePath] = uigetfile([saveLoc, '/Numerical/MATLAB/planarSprayPOD/*.mat'], ...
                                             'Select POD Data');

            switch format

                case 'A'

                    if contains(filePath, '/base')
                        disp(['Loading ''', fileName, '''...']);
                        
                        campaignID = load([filePath, fileName], 'campaignID').campaignID;
                        caseID = load([filePath, fileName], 'caseID').caseID;
                        dataID = load([filePath, fileName], 'dataID').dataID;
                        PODdata = load([filePath, fileName], 'PODdata').PODdata;
                        sampleInt = load([filePath, fileName], 'sampleInt').sampleInt;
                        timePrecision = load([filePath, fileName], 'timePrecision').timePrecision;
                        dLims = load([filePath, fileName], 'dLims').dLims;
                        
                        disp('    Success');
                        
                        valid = true;
                    else
                        disp('WARNING: Invalid File Selection');
                        clear fileName filePath;
                    end

                case 'B'

                    if contains(filePath, '/X_')
                        disp(['Loading ''', fileName, '''...']);
                        
                        campaignID = load([filePath, fileName], 'campaignID').campaignID;
                        caseID = load([filePath, fileName], 'caseID').caseID;
                        planeID = load([filePath, fileName], 'planeID').planeID;
                        dataID = load([filePath, fileName], 'dataID').dataID;
                        PODdata = load([filePath, fileName], 'PODdata').PODdata;
                        sampleInt = load([filePath, fileName], 'sampleInt').sampleInt;
                        timePrecision = load([filePath, fileName], 'timePrecision').timePrecision;
                        dLims = load([filePath, fileName], 'dLims').dLims;
                        
                        disp('    Success');

                        valid = true;
                    else
                        disp('WARNING: Invalid File Selection');
                        clear fileName filePath;
                    end

            end
        end
        clear valid;
        
    case 'C'
        
        valid = false;
        while ~valid
            disp(' ');
            
            [fileName, filePath] = uigetfile([saveLoc, '/Experimental/MATLAB/planarSprayPOD/*.mat'], ...
                                             'Select POD Data');
            
            if contains(filePath, 'Hz')
                disp(['Loading ''', fileName, '''...']);
                
                campaignID = load([filePath, fileName], 'campaignID').campaignID;
                caseID = load([filePath, fileName], 'caseID').caseID;
                planeID = load([filePath, fileName], 'planeID').planeID;
                dataID = load([filePath, fileName], 'dataID').dataID;
                PODdata = load([filePath, fileName], 'PODdata').PODdata;
                cellSize = load([filePath, fileName], 'cellSize').cellSize;
                sampleInt = load([filePath, fileName], 'sampleInt').sampleInt;
                timePrecision = load([filePath, fileName], 'timePrecision').timePrecision;
                
                disp('    Success');
                
                valid = true;
            else
                disp('WARNING: Invalid File Selection');
                clear fileName filePath;
            end
            
        end
        clear valid;
        
end

field = PODdata.POD.field;

Ns = height(PODdata.positionGrid);
Nt = height(PODdata.time);

modes90 = find(cumsum(PODdata.POD.modeEnergy) > 90, 1);
modes99 = find(cumsum(PODdata.POD.modeEnergy) > 99, 1);

if strcmp(campaignID, 'Far_Field_Soiling_07_22')
    refValue = 0.0052166;
end

if ~exist('refValue', 'var')
    refValue = 1;
end

disp(' ');
disp(' ');


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry(geoLoc);

disp(' ');
disp(' ');


%% Select Reconstruction Options

disp('Reconstruction Options');
disp('-----------------------');

disp(' ');

disp(['Selected Decomposition Contains ', num2str(Nt), ' Modes:']);
disp(['    The First ', num2str(modes90), ' Modes Contain Approximately 90% of the Total Captured Energy Content']);
disp(['    The First ', num2str(modes99), ' Modes Contain Approximately 99% of the Total Captured Energy Content']);


valid = false;
while ~valid
    disp(' ');
    
    nModes = inputModes(Nt);

    if nModes == -1
        continue;
    else
        nModes = sort(nModes);
    end

    valid = true;
end
clear valid;

disp(' ');
disp(' ');


%% Perform Field Reconstruction

disp('Field Reconstruction');
disp('---------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

disp(' ');

disp('    Initialising...');

% Specify Map Boundaries
switch format
    
    case 'A'
        % Identify Model Base
        parts = fieldnames(geometry);
        for i = 1:height(parts)
            
            if max(geometry.(parts{i}).vertices(:,1)) == xDims(2)
                parts = parts{i};
                break;
            end
            
            if i == height(parts)
                error('Mismatch Between ''xDims'' and Geometry Bounding Box')
            end
            
        end
        clear i;
    
        geoPoints = geometry.(parts).vertices;
        basePoints = geoPoints((geoPoints(:,1) == xDims(2)),:);
        
        mapPerim = boundary(basePoints(:,2), basePoints(:,3), 1);
        mapPerim = basePoints(mapPerim,:);
        basePoly = polyshape(mapPerim(:,2), mapPerim(:,3), 'keepCollinearPoints', true);
        basePoly = polybuffer(basePoly, -0.005, 'jointType', 'square');
        mapPerim = ones(height(basePoly.Vertices),3) * mapPerim(1,1);
        mapPerim(:,[2,3]) = basePoly.Vertices(:,[1,2]);

        if ~all(mapPerim(1,:) == mapPerim(end,:))
            mapPerim = [mapPerim; mapPerim(1,:)]; % Close Boundary
        end
        
        xLimsData = xDims(2) + 1e-3; % Offset Particles From Base for Better Visibility
        yLimsData = [min(mapPerim(:,2)); max(mapPerim(:,2))];
        zLimsData = [min(mapPerim(:,3)); max(mapPerim(:,3))];
        
    case {'B', 'C'}
        mapPerim = [];

        xLimsData = PODdata.positionGrid(1,1);
        yLimsData = [min(PODdata.positionGrid(:,2)); max(PODdata.positionGrid(:,2))];
        zLimsData = [min(PODdata.positionGrid(:,3)); max(PODdata.positionGrid(:,3))];
        
end

% Initialise Reconstruction Variables
reconData = PODdata;
reconData.(field) = []; reconData.(field).true = PODdata.(field);

if (any(strcmp(format, {'A', 'B'})) && any(strcmp(field, 'mass'))) || strcmp(format, 'C')
    reconData.CoM = []; reconData.CoM.true = PODdata.CoM;
end

clear PODdata;

for i = 1:Nt
    reconData.(field).recon.inst{i} = reconData.(field).true.mean;
end
clear i;

disp(' ');

% Perform Field Reconstruction
if width(nModes) == 1
    reconData = reconstructPOD(reconData, field, nModes, Ns, Nt, 'scalar', true);
else
    reconData = reconstructPOD(reconData, field, nModes, Ns, Nt, 'scalar', false);
end

disp(' ');

% Calculate Calculate Reconstructed Time Average
disp('    Calculating Reconstructed Time Average');

reconData.(field).recon.mean = zeros([Ns,1]);

for i = 1:Nt
    reconData.(field).recon.mean = reconData.(field).recon.mean + reconData.(field).recon.inst{i};
end
clear i;

reconData.(field).recon.mean = reconData.(field).recon.mean / Nt;

disp(' ');

% Calculate Instantaneous Spray Density Fluctuations
disp('    Calculating Reconstructed Instantaneous Fluctuations');

reconData.(field).recon.prime = reconData.(field).recon.inst;

for i = 1:Nt
    reconData.(field).recon.prime{i} = reconData.(field).recon.inst{i} - reconData.(field).recon.mean;
end
clear i;

disp(' ');

% Calculate RMS of Reconstructed Field
disp('    Calculating RMS of Reconstructed Field');

reconData.(field).recon.RMS = zeros([Ns,1]);

for i = 1:Nt
    reconData.(field).recon.RMS  = reconData.(field).recon.RMS + reconData.(field).recon.prime{i}.^2;
end
clear i;

reconData.(field).recon.RMS = sqrt((1 / Nt) * reconData.(field).recon.RMS);

if (any(strcmp(format, {'A', 'B'})) && any(strcmp(field, 'mass'))) || strcmp(format, 'C')
    disp(' ');

    disp('    Calculating Reconstructed Centre of Mass...');

    % Initialise Progress Bar
    wB = waitbar(0, 'Calculating Reconstructed Centre of Mass', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    
    % Calculate Reconstructed CoM
    reconData.CoM.recon.inst = cell(Nt,1); reconData.CoM.recon.inst(:) = {zeros([1,3])};
    
    for i = 1:Nt
        reconData.CoM.recon.inst{i}(1) = reconData.positionGrid(1,1);
        reconData.CoM.recon.inst{i}(2) = sum(reconData.(field).recon.inst{i} .* reconData.positionGrid(:,2)) / sum(reconData.(field).recon.inst{i});
        reconData.CoM.recon.inst{i}(3) = sum(reconData.(field).recon.inst{i} .* reconData.positionGrid(:,3)) / sum(reconData.(field).recon.inst{i});
        
        % Update Waitbar
        waitbar((i / Nt), wB);
    end
    clear i;
    
    delete(wB);
    
    disp(' ');

    % Calculate Time-Averaged Centre of Spray
    disp('    Calculating Time-Averaged Centre of Spray');

    reconData.CoM.recon.mean = zeros([1,3]);

    reconData.CoM.recon.mean(1) = reconData.positionGrid(1,1);
    reconData.CoM.recon.mean(2) = sum(reconData.(field).recon.mean .* reconData.positionGrid(:,2)) / sum(reconData.(field).recon.mean);
    reconData.CoM.recon.mean(3) = sum(reconData.(field).recon.mean .* reconData.positionGrid(:,3)) / sum(reconData.(field).recon.mean);
end

%%%%

executionTime = toc;

disp(' ');

disp(['    Run Time: ', num2str(executionTime), 's']);

disp(' ');

disp('  SUCCESS  ');
disp('***********');

disp(' ');
disp(' ');


%% Select Reconstruction Presentation Options

disp('Reconstruction Presentation Options');
disp('------------------------------------');

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Plot Average Reconstructed Field? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotMean = false;
        
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotMean = true;
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Plot RMS of Reconstructed Field? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotRMS = false;
        
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotRMS = true;
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Plot Instantaneous Reconstructed Fields? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotInst = false;
        
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotInst = true;
        
        startFrame = inputFrames(Nt, 'Start');
        
        if startFrame == -1
            continue;
        end
        
        endFrame = inputFrames(Nt, 'End');
        
        if endFrame == -1
            continue;
        elseif endFrame < startFrame
            disp('        WARNING: Invalid Time Format (''endFrame'' Precedes ''startFrame'')');
            
            continue;
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

if plotMean || plotRMS || plotInst
    
    % Normalise Coordinate System
    if normDims
        disp(' ');

        disp('Normalising Spatial Dimensions...');

        wB = waitbar(0, 'Normalising Spatial Dimensions', 'name', 'Progress');
        wB.Children.Title.Interpreter = 'none';

        parts = fieldnames(geometry);
        for i = 1:height(parts)
            geometry.(parts{i}).vertices = geometry.(parts{i}).vertices / normLength;
        end
        clear i parts;

        xDims = xDims / normLength;
        yDims = yDims / normLength;
        zDims = zDims / normLength;

        mapPerim = mapPerim / normLength;

        xLimsData = xLimsData / normLength;
        yLimsData = yLimsData / normLength;
        zLimsData = zLimsData / normLength;

        cellSize.target = cellSize.target / normLength;
        cellSize.y = cellSize.y / normLength;
        cellSize.z = cellSize.z / normLength;
        cellSize.area = cellSize.area / (normLength^2);

        reconData.positionGrid = reconData.positionGrid / normLength;

        reconData.CoM.recon.mean = reconData.CoM.recon.mean / normLength;

        for i = 1:Nt
            reconData.CoM.recon.inst{i} = reconData.CoM.recon.inst{i} / normLength;

            % Update Waitbar
            waitbar((i / Nt), wB);
        end
        clear i;

        delete(wB);
    end
    
    % Normalise Spray Density
    if normDensity
        disp(' ');
        
        disp('Normalising Spray Density...');
        
        wB = waitbar(0, 'Normalising Spray Density', 'name', 'Progress');
        wB.Children.Title.Interpreter = 'none';
        
        reconData.density.recon.mean = reconData.density.recon.mean / refValue;
        mapData.density.recon.RMS = reconData.density.recon.RMS / refValue;
        
        for i = 1:Nt
            reconData.density.recon.inst{i} = reconData.density.recon.inst{i} / refValue;
            reconData.density.recon.prime{i} = reconData.density.recon.prime{i} / refValue;
            
            % Update Waitbar
            waitbar((i / Nt), wB);
        end
        clear i;
        
        delete(wB);
    end
    
    disp(' ');
    
    % Select Variable(s) of Interest
    disp('Select Variable(s) To Plot...');

    dataVars = fieldnames(reconData);
    nonFieldVars = {'positionGrid'; 'time'; 'CoM'; 'POD'};
    fieldVars = setdiff(dataVars, nonFieldVars);
    clear dataVars nonFieldVars;

    valid = false;
    while ~valid
        [index, valid] = listdlg('listSize', [300, 300], ...
                                 'selectionMode', 'multiple', ...
                                 'name', 'Select Variable(s) to Plot', ...
                                 'listString', fieldVars);

        if ~valid
            disp('    WARNING: No Valid Variable Selected');
        end

    end
    clear valid;

    plotVars = fieldVars(index); clear fieldVars;
    
    disp(['    Plotting ', num2str(height(plotVars)), ' Variable(s) of Interest']);
end

disp(' ');
disp(' ');


%% Present Reconstruction

disp('Reconstruction Presentation');
disp('----------------------------');

disp(' ');

if plotMean || plotRMS || plotInst
    orientation = 'YZ';
    positionData = reconData.positionGrid;
    
    if strcmp(campaignID, 'Windsor_fullScale')
        spatialRes = 2e-3;
    elseif strcmp(campaignID, 'Windsor_Upstream_2023')
        spatialRes = 0.5e-3;
    else
        spatialRes = 0.5e-3;
    end
    
    if normDims
        spatialRes = spatialRes / normLength;
    end
    
    % Offset Particles From Surface to Improve Visibility
    switch format
        
        case 'A'
            xLimsData = xLimsData + spatialRes;
            positionData(:,1) = xLimsData;
            
    end
    
    nPlanes = 1;
    planeNo = 1;
    cMap = flipud(viridis(32));
    refPoint = [];
    
    switch format

        case 'A'
            xLimsPlot = [0.3; 4.6257662];
            yLimsPlot = [-0.25; 0.25];
            zLimsPlot = [0; 0.4];
            
        case {'B', 'C'}
            xLimsPlot = [0.3; 4.6257662];
            yLimsPlot = [-0.5; 0.5];
            zLimsPlot = [0; 0.5];
            
    end
    
    if ~normDims
        xLimsPlot = xLimsPlot * normLength;
        yLimsPlot = yLimsPlot * normLength;
        zLimsPlot = zLimsPlot * normLength;
    end
    
end

if plotMean
    disp('    Presenting Average Reconstructed  Field...');

    scalarData = reconData.(field).recon.mean;
    figName = ['Recon_Average_', caseID];
    figTitle = '{ }'; % Leave Blank ('-') for Formatting Purposes
    
    switch format
        
        case {'A', 'B'}
            contourlines = [];
            cLims = [0; max(cellfun(@max, reconData.(field).recon.inst))];
            
        case 'C'
            contourlines = [0.02; 0.02];
            cLims = [0; 1.05];
    
    end

    [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                           xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                           planeNo, fig, figName, cMap, geometry, contourlines, ...
                                           refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                           zLimsPlot, normDims, figSave);
    
    disp(' ');
end

if plotRMS
    disp('    Presenting RMS of Reconstructed Field...');
    
    scalarData = reconData.(field).recon.RMS;
    figName = ['Recon_RMS_', caseID];
    contourlines = [];
    figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
    
    switch format
        
        case {'A', 'B'}
            contourlines = [];
            cLims = [0; max(cellfun(@max, reconData.(field).recon.inst))];
            
        case 'C'
            contourlines = [];
            cLims = [0; 5];
    
    end

    [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                           xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                           planeNo, fig, figName, cMap, geometry, contourlines, ...
                                           refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                           zLimsPlot, normDims, figSave);
    
    disp(' ');
end

if plotInst
    disp('    Presenting Instantaneous Reconstructed Field...');
    
    switch format
        
        case {'A', 'B'}
            cLims = [0; max(cellfun(@max, reconData.(field).recon.inst))];
            contourlines = [];
            
        case 'C'
            cLims = [0; 3.6];
            contourlines = [];
            
    end

    figHold = fig;
    
    for i = startFrame:endFrame
        
        if i ~= startFrame
            clf(fig);
            fig = figHold;
        end
        
        scalarData = reconData.(field).recon.inst{i};
        figTime = num2str(reconData.time(i), ['%.', num2str(timePrecision), 'f']);
        
        switch format
            
            case 'A'
                figName = ['Base_', field, '_Recon_T', erase(figTime, '.'), '_', caseID];
            
            case 'B'
                figName = [planeID, '_', field, '_Recon_T', erase(figTime, '.'), '_', caseID];
                
            case 'C'
                figName = ['Recon_T', erase(figTime, '.'), '_', caseID];
        
        end
        
        figTitle = ['{', figTime, ' \it{s}}'];
        
        [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                               xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                               planeNo, fig, figName, cMap, geometry, contourlines, ...
                                               refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                               zLimsPlot, normDims, figSave);
    end
    clear i;
    
    disp(' ');
end

if ~plotMean && ~plotRMS && ~plotInst
    disp('    Skipping Map Presentation');
    
    disp(' ');
end

disp(' ');


%% Save Reconstruction Data

disp('Data Save Options');
disp('------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Save Data for Future Use? [y/n]: ', 's');
    
    if selection == 'n' | selection == 'N' %#ok<OR2>
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        
        switch format
            
            case 'A'
                
                if ~exist([saveLoc, '/Numerical/MATLAB/planarSprayReconstruction/', campaignID, '/', caseID, '/base/', field], 'dir')
                    mkdir([saveLoc, '/Numerical/MATLAB/planarSprayReconstruction/', campaignID, '/', caseID, '/base/', field]);
                end
                
            case 'B'
                
                if ~exist([saveLoc, '/Numerical/MATLAB/planarSprayReconstruction/', campaignID, '/', caseID, '/', planeID, '/', field], 'dir')
                    mkdir([saveLoc, '/Numerical/MATLAB/planarSprayReconstruction/', campaignID, '/', caseID, '/', planeID, '/', field]);
                end
                
            case 'C'
                
                if ~exist([saveLoc, '/Experimental/MATLAB/planarSprayReconstruction/', campaignID, '/', caseID, '/', field], 'dir')
                    mkdir([saveLoc, '/Experimental/MATLAB/planarSprayReconstruction/', campaignID, '/', caseID, '/', field]);
                end
                
        end
        
        switch format
            
            case 'A'
                disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/planarSprayReconstruction/', campaignID, '/', caseID, '/base/', field, '/', dataID, '_', mat2str(nModes), '.mat']);
                save([saveLoc, '/Numerical/MATLAB/planarSprayReconstruction/', campaignID, '/', caseID, '/base/', field, '/', dataID, '_', mat2str(nModes), '.mat'], ...
                     'campaignID', 'caseID', 'dataID', 'reconData', 'nModes', 'cellSize', 'sampleInt', 'dLims', 'timePrecision', 'normDims', '-v7.3', '-noCompression');
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/planarSprayReconstruction/', campaignID, '/', caseID, '/', planeID, '/', field, '/', dataID, '_', mat2str(nModes), '.mat']);
                save([saveLoc, '/Numerical/MATLAB/planarSprayReconstruction/', campaignID, '/', caseID, '/', planeID, '/', field, '/', dataID, '_', mat2str(nModes), '.mat'], ...
                     'campaignID', 'caseID', 'planeID', 'dataID', 'reconData', 'nModes', 'cellSize', 'sampleInt', 'dLims', 'timePrecision', 'normDims', '-v7.3', '-noCompression');
                disp('        Success');
                
            case 'C'
                disp(['    Saving to: ', saveLoc, '/Experimental/MATLAB/planarSprayReconstruction/', campaignID, '/', caseID, '/', field, '/', dataID, '_', mat2str(nModes), '.mat']);
                save([saveLoc, '/Experimental/MATLAB/planarSprayPOD/', campaignID, '/', caseID, '/', field, '/', dataID, '_', mat2str(nModes), '.mat'], ...
                     'campaignID', 'caseID', 'planeID', 'dataID', 'reconData', 'nModes', 'cellSize', 'sampleInt', 'timePrecision', 'normDims', '-v7.3', '-noCompression');
                disp('        Success');
        
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;


%% Local Functions

function modes = inputModes(Nt)

    modes = str2num(input('    Input Select the Desired Modes To Reconstruct [Row Vector Form]: ', 's')); %#ok<ST2NM>
    
    if isempty(modes) || any(isnan(modes)) || ~isrow(modes) > 1 || any(modes <= 0) || any(modes > Nt)
        disp('        WARNING: Invalid Entry');
        
        modes = -1;
    end

end


function frameNo = inputFrames(Nt, type)

    frameNo = str2double(input(['    Input Desired ', type, ' Frame [1-', num2str(Nt), ']: '], 's'));
    
    if isnan(frameNo) || frameNo < 1 || frameNo > Nt
        disp('        WARNING: Invalid Entry');
        
        frameNo = -1;
    end

end