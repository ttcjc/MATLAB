%% Preamble

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

if exist('/mnt/Processing/Data', 'dir')
    saveLocation = '/mnt/Processing/Data';
else
    saveLocation = '~/Data';
end

nProc = maxNumCompThreads - 2; % Number of Processors Used for Process-Based Parallelisation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures


%% Planar Spray POD Field Reconstructor v1.0

figSave = false; % Save .fig File(s)

disp('=========================================');
disp('Planar Spray POD Field Reconstructor v1.0');
disp('=========================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit (Functionality Separated From 'planarSprayPOD' To Improve Efficiency)


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
            
            [fileName, filePath] = uigetfile([saveLocation, '/Numerical/MATLAB/planarSprayPOD/*.mat'], ...
                                             'Select POD Data');

            switch format

                case 'A'

                    if contains(filePath, '/base')
                        disp(['Loading ''', fileName, '''...']);
                        
                        caseID = load([filePath, fileName], 'caseID').caseID;
                        dataID = load([filePath, fileName], 'dataID').dataID;
                        PODdata = load([filePath, fileName], 'PODdata').PODdata;
                        sampleInt = load([filePath, fileName], 'sampleInt').sampleInt;
                        timePrecision = load([filePath, fileName], 'timePrecision').timePrecision;
                        dLims = load([filePath, fileName], 'dLims').dLims;
                        normDims = load([filePath, fileName], 'normDims').normDims;
                        
                        disp('    Success');
                        
                        valid = true;
                    else
                        disp('WARNING: Invalid File Selection');
                        clear fileName filePath;
                    end

                case 'B'

                    if contains(filePath, '/X_')
                        disp(['Loading ''', fileName, '''...']);
                        
                        caseID = load([filePath, fileName], 'caseID').caseID;
                        planeID = load([filePath, fileName], 'planeID').planeID;
                        dataID = load([filePath, fileName], 'dataID').dataID;
                        PODdata = load([filePath, fileName], 'PODdata').PODdata;
                        sampleInt = load([filePath, fileName], 'sampleInt').sampleInt;
                        timePrecision = load([filePath, fileName], 'timePrecision').timePrecision;
                        dLims = load([filePath, fileName], 'dLims').dLims;
                        normDims = load([filePath, fileName], 'normDims').normDims;
                        
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
            
            [fileName, filePath] = uigetfile([saveLocation, '/Experimental/MATLAB/planarSprayPOD/*.mat'], ...
                                             'Select POD Data');
            
            if contains(filePath, 'Hz')
                disp(['Loading ''', fileName, '''...']);
                
                caseID = load([filePath, fileName], 'caseID').caseID;
                planeID = load([filePath, fileName], 'planeID').planeID;
                dataID = load([filePath, fileName], 'dataID').dataID;
                PODdata = load([filePath, fileName], 'PODdata').PODdata;
                cellSize = load([filePath, fileName], 'cellSize').cellSize;
                sampleInt = load([filePath, fileName], 'sampleInt').sampleInt;
                timePrecision = load([filePath, fileName], 'timePrecision').timePrecision;
                normDims = load([filePath, fileName], 'normDims').normDims;
                
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

disp(' ');
disp(' ');


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normDims, normLength] = selectGeometry(normDims);

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


disp(' ');
disp(' ');


%% Present Reconstruction

disp('Reconstruction Presentation');
disp('----------------------------');

disp(' ');

if plotMean || plotRMS || plotInst
    orientation = 'YZ';
    
    switch format

        case 'A'
            
            if contains(caseID, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 1.43075];
                yLimsPlot = [-0.2445; 0.2445];
                zLimsPlot = [0; 0.389];
            end

        case 'B'
            
            if contains(caseID, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 2.73575];
                yLimsPlot = [-0.522; 0.522];
                zLimsPlot = [0; 0.6264];
            end
            
        case 'C'
            xLimsPlot = [0.31875; 2.73575];
            yLimsPlot = [-0.522; 0.522];
            zLimsPlot = [0; 0.522];
            
    end

    if normDims
        xLimsPlot = round((xLimsPlot / normLength), spacePrecision);
        yLimsPlot = round((yLimsPlot / normLength), spacePrecision);
        zLimsPlot = round((zLimsPlot / normLength), spacePrecision);
    end
    
    switch format
        
        case {'A', 'B'}
            
            if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
                spatialRes = 0.5e-3;
            else
                spatialRes = 2e-3;
            end
            
        case 'C'
            spatialRes = 0.5e-3;
            
    end

    positionData = reconData.positionGrid;
    nPlanes = 1;
    planeNo = 1;
    cMap = flipud(viridis(32));
    refPoint = [];
    
    % Remove Geometry From Empty Tunnel
    if contains(caseID, 'ET')
        geometry = [];
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
    cLims = 'auto';

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
                
                if ~exist([saveLocation, '/Numerical/MATLAB/planarSprayReconstruction/', caseID, '/base/', field], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/planarSprayReconstruction/', caseID, '/base/', field]);
                end
                
            case 'B'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/planarSprayReconstruction/', caseID, '/', planeID, '/', field], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/planarSprayReconstruction/', caseID, '/', planeID, '/', field]);
                end
                
            case 'C'
                
                if ~exist([saveLocation, '/Experimental/MATLAB/planarSprayReconstruction/', caseID, '/', field], 'dir')
                    mkdir([saveLocation, '/Experimental/MATLAB/planarSprayReconstruction/', caseID, '/', field]);
                end
                
        end
        
        switch format
            
            case 'A'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/planarSprayReconstruction/', caseID, '/base/', field, '/', dataID, '_', mat2str(nModes), '.mat']);
                save([saveLocation, '/Numerical/MATLAB/planarSprayReconstruction/', caseID, '/base/', field, '/', dataID, '_', mat2str(nModes), '.mat'], ...
                     'caseID', 'dataID', 'reconData', 'nModes', 'cellSize', 'sampleInt', 'dLims', 'timePrecision', 'normDims', '-v7.3', '-noCompression');
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/planarSprayReconstruction/', caseID, '/', planeID, '/', field, '/', dataID, '_', mat2str(nModes), '.mat']);
                save([saveLocation, '/Numerical/MATLAB/planarSprayReconstruction/', caseID, '/', planeID, '/', field, '/', dataID, '_', mat2str(nModes), '.mat'], ...
                     'caseID', 'planeID', 'dataID', 'reconData', 'nModes', 'cellSize', 'sampleInt', 'dLims', 'timePrecision', 'normDims', '-v7.3', '-noCompression');
                disp('        Success');
                
            case 'C'
                disp(['    Saving to: ', saveLocation, '/Experimental/MATLAB/planarSprayReconstruction/', caseID, '/', field, '/', dataID, '_', mat2str(nModes), '.mat']);
                save([saveLocation, '/Experimental/MATLAB/planarSprayPOD/', caseID, '/', field, '/', dataID, '_', mat2str(nModes), '.mat'], ...
                     'caseID', 'planeID', 'dataID', 'reconData', 'nModes', 'cellSize', 'sampleInt', 'timePrecision', 'normDims', '-v7.3', '-noCompression');
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