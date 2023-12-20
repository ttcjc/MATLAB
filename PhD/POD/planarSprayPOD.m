%% Planar Spray POD Calculator v3.2
% ----
% Perform Proper Orthogonal Decomposition on Previously Processed Planar Spray Maps
% (Collected Using 'planarSprayMapper' & 'expPlanarSprayMapper')


%% Preamble

run preamble;

normDims = true; % Normalise Spatial Dimensions in Plots

flipMode = true; % Present Both Orientations of Mode(s)

figSave = false; % Save .fig File(s)

disp('================================');
disp('Planar Spray POD Calculator v3.2');
disp('================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit (Base Contamination and Far-Field Extraction Plane)
% v1.1 - Rename and Restructure to Account for Changes to 'mapData' Output
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats
% v2.1 - Moved Calculation of Instantaneous Variables To Pre-Processing Scripts
% v3.0 - Offloaded Reconstruction To Improve Efficiency
% v3.1 - Minor Update to Shift Preamble Into Separate Script
% v3.2 - Update To Correct Inconsistent Normalisation Throughout Repository


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

disp('Spray Map Acquisition');
disp('----------------------');

switch format
    
    case {'A', 'B'}
        
        valid = false;
        while ~valid
            disp(' ');
            
            [fileName, filePath] = uigetfile([saveLoc, '/Numerical/MATLAB/planarSprayMap/*.mat'], ...
                                             'Select Map Data');

            switch format

                case 'A'

                    if contains(filePath, '/base')
                        disp(['Loading ''', fileName, '''...']);
                        
                        campaignID = load([filePath, fileName], 'campaignID').campaignID;
                        caseID = load([filePath, fileName], 'caseID').caseID;
                        dataID = load([filePath, fileName], 'dataID').dataID;
                        mapData = load([filePath, fileName], 'mapData').mapData;
                        cellSize = load([filePath, fileName], 'cellSize').cellSize;
                        sampleInt = load([filePath, fileName], 'sampleInt').sampleInt;
                        timePrecision = load([filePath, fileName], 'timePrecision').timePrecision;
                        dLims = load([filePath, fileName], 'dLims').dLims;
                        dataFormat = load([filePath, fileName], 'dataFormat').dataFormat;
                        
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
                        mapData = load([filePath, fileName], 'mapData').mapData;
                        cellSize = load([filePath, fileName], 'cellSize').cellSize;
                        sampleInt = load([filePath, fileName], 'sampleInt').sampleInt;
                        timePrecision = load([filePath, fileName], 'timePrecision').timePrecision;
                        dLims = load([filePath, fileName], 'dLims').dLims;
                        dataFormat = load([filePath, fileName], 'dataFormat').dataFormat;
                        
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
            
            [fileName, filePath] = uigetfile([saveLoc, '/Experimental/MATLAB/planarSprayMap/*.mat'], ...
                                             'Select Map Data');
            
            if contains(filePath, 'Hz')
                disp(['Loading ''', fileName, '''...']);
                
                campaignID = load([filePath, fileName], 'campaignID').campaignID;
                caseID = load([filePath, fileName], 'caseID').caseID;
                planeID = load([filePath, fileName], 'planeID').planeID;
                dataID = load([filePath, fileName], 'dataID').dataID;
                mapData = load([filePath, fileName], 'mapData').mapData;
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

disp(' ');
disp(' ');


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry;

disp(' ');
disp(' ');


%% Select POD Options

disp('POD Options');
disp('------------');

disp(' ');

% Select Variable of Interest
disp('Select Variable for Decomposition...');

mapDataVars = fieldnames(mapData);
nonFieldVars = {'positionGrid'; 'time'; 'CoM'};
fieldVars = setdiff(mapDataVars, nonFieldVars);
clear mapDataVars nonFieldVars;

valid = false;
while ~valid
    [index, valid] = listdlg('listSize', [300, 300], ...
                             'selectionMode', 'single', ...
                             'name', 'Select Variable for Decomposition', ...
                             'listString', fieldVars);
    
    if valid && ~isfield(mapData.(fieldVars{index}), 'prime')
        disp('    WARNING: Selected Variable Is Not a Valid Field for Decomposition');
        
        valid = false;        
        continue;
    end

    if ~valid
        disp('    WARNING: No Mapping Variable Selected');
    end
    
end
clear valid;

field = fieldVars{index}; clear fieldVars;

disp(['    Variable of Interest: ', field]);

disp(' ');
disp(' ');


%% Perform Planar POD (Snapshot Method)

disp('Planar Proper Orthogonal Decomposition');
disp('---------------------------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

%%%%

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
        
        mapPerim = boundary(basePoints(:,2), basePoints(:,3), 0.95);
        mapPerim = basePoints(mapPerim,:);
        basePoly = polyshape(mapPerim(:,2), mapPerim(:,3), 'keepCollinearPoints', true);

        if strcmp(campaignID, 'Windsor_fullScale')
            basePoly = polybuffer(basePoly, -16e-3, 'jointType', 'square');
        elseif strcmp(campaignID, 'Windsor_Upstream_2023')
            basePoly = polybuffer(basePoly, -4e-3, 'jointType', 'square');
        else
            basePoly = polybuffer(basePoly, -4e-3, 'jointType', 'square');
        end
        
        mapPerim = ones(height(basePoly.Vertices),3) * mapPerim(1,1);
        mapPerim(:,[2,3]) = basePoly.Vertices(:,[1,2]);

        if ~all(mapPerim(1,:) == mapPerim(end,:))
            mapPerim = [mapPerim; mapPerim(1,:)]; % Close Boundary
        end
        
        xLimsData = xDims(2);
        yLimsData = [min(mapPerim(:,2)); max(mapPerim(:,2))];
        zLimsData = [min(mapPerim(:,3)); max(mapPerim(:,3))];
        
    case {'B', 'C'}
        mapPerim = [];

        xLimsData = mapData.positionGrid(1,1);
        yLimsData = [min(mapData.positionGrid(:,2)); max(mapData.positionGrid(:,2))];
        zLimsData = [min(mapData.positionGrid(:,3)); max(mapData.positionGrid(:,3))];
        
end

% Initialise POD Variables
PODdata = mapData; clear mapData;

disp(' ');

% Perform Planar Snapshot POD
switch format

    case 'A'
        [fig, PODdata, modesEnergetic, modes80percent, Ns, Nt] = performPOD(fig, PODdata, field, ...
                                                                            'scalar', 'Base', figSave);

    case 'B'
        [fig, PODdata, modesEnergetic, modes80percent, Ns, Nt] = performPOD(fig, PODdata, field, ...
                                                                            'scalar', planeID, figSave);
    
    case 'C'
        [fig, PODdata, modesEnergetic, modes80percent, Ns, Nt] = performPOD(fig, PODdata, field, ...
                                                                            'scalar', planeID, figSave);

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


%% Save POD Data

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
                
                if ~exist([saveLoc, '/Numerical/MATLAB/planarSprayPOD/', campaignID, '/', caseID, '/base/', field], 'dir')
                    mkdir([saveLoc, '/Numerical/MATLAB/planarSprayPOD/', campaignID, '/', caseID, '/base/', field]);
                end
                
            case 'B'
                
                if ~exist([saveLoc, '/Numerical/MATLAB/planarSprayPOD/', campaignID, '/', caseID, '/', planeID, '/', field], 'dir')
                    mkdir([saveLoc, '/Numerical/MATLAB/planarSprayPOD/', campaignID, '/', caseID, '/', planeID, '/', field]);
                end
                
            case 'C'
                
                if ~exist([saveLoc, '/Experimental/MATLAB/planarSprayPOD/', campaignID, '/', caseID, '/', field], 'dir')
                    mkdir([saveLoc, '/Experimental/MATLAB/planarSprayPOD/', campaignID, '/', caseID, '/', field]);
                end
                
        end
        
        switch format
            
            case 'A'
                disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/planarSprayPOD/', campaignID, '/', caseID, '/base/', field, '/', dataID, '.mat']);
                save([saveLoc, '/Numerical/MATLAB/planarSprayPOD/', campaignID, '/', caseID, '/base/', field, '/', dataID, '.mat'], ...
                      'campaignID', 'caseID', 'dataID', 'PODdata', 'cellSize', 'sampleInt', 'timePrecision', 'dataFormat', '-v7.3', '-noCompression');
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/planarSprayPOD/', campaignID, '/', caseID, '/', planeID, '/', field, '/', dataID, '.mat']);
                save([saveLoc, '/Numerical/MATLAB/planarSprayPOD/', campaignID, '/', caseID, '/', planeID, '/', field, '/', dataID, '.mat'], ...
                      'campaignID', 'caseID', 'planeID', 'dataID', 'PODdata', 'cellSize', 'sampleInt', 'timePrecision', 'dataFormat', '-v7.3', '-noCompression');
                disp('        Success');
            
            case 'C'
                disp(['    Saving to: ', saveLoc, '/Experimental/MATLAB/planarSprayPOD/', campaignID, '/', caseID, '/', field, '/', dataID, '.mat']);
                save([saveLoc, '/Experimental/MATLAB/planarSprayPOD/', campaignID, '/', caseID, '/', field, '/', dataID, '.mat'], ...
                      'campaignID', 'caseID', 'planeID', 'dataID', 'PODdata', 'cellSize', 'sampleInt', 'timePrecision', '-v7.3', '-noCompression');
                disp('        Success');
        
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Select Mode Presentation Options

disp('Mode Presentation Options');
disp('--------------------------');

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Plot Individual Modes? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotModes = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotModes = true;
        nModes = inputModes(Nt);
        
        if nModes == -1
            continue;
        else
            nModes = sort(nModes);
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

if plotModes
    
    % Normalise Coordinate System
    if normDims
        disp(' ');

        disp('    Normalising Spatial Dimensions...');

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
        
        switch format
            
            case 'C'
                cellSize.DaVis = cellSize.DaVis / normLength;
            
        end

        cellSize.target = cellSize.target / normLength;
        cellSize.y = cellSize.y / normLength;
        cellSize.z = cellSize.z / normLength;
        cellSize.area = cellSize.area / (normLength^2);

        PODdata.positionGrid = PODdata.positionGrid / normLength;

        PODdata.CoM.mean = PODdata.CoM.mean / normLength;

        for i = 1:Nt
            PODdata.CoM.inst{i} = PODdata.CoM.inst{i} / normLength;

            % Update Waitbar
            waitbar((i / Nt), wB);
        end
        clear i;

        delete(wB);
    end
    
end

disp(' ');
disp(' ');


%% Present POD Modes

disp('Mode Presentation');
disp('------------------');

disp(' ');

if plotModes
    orientation = 'YZ';
    positionData = PODdata.positionGrid;
    
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
    cMap = cool2warm(32);
    contourlines = [];
    refPoint = [];
    cLims = [-1; 1];
    
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
    
    for i = nModes
        disp(['    Presenting Mode #', num2str(i), '...']);

        scalarData = rescale(PODdata.POD.phi(:,i), -1, 1);

        switch format

            case 'A'
                figName = ['Base_POD_', field, '_M', num2str(i)];

            case 'B'
                figName = [planeID, '_POD_', field, '_M', num2str(i)];
            
            case 'C'
                figName = ['POD_M', num2str(i), '_', caseID];

        end

        figTitle = [num2str(round(PODdata.POD.modeEnergy(i), 2), '%.2f'), '%'];

        [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                               xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                               planeNo, fig, figName, cMap, geometry, contourlines, ...
                                               refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                               zLimsPlot, normDims, figSave);
        
        if flipMode
            [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                                   xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                                   planeNo, fig, [figName, '_Flip'], flipud(cMap), ...
                                                   geometry, contourlines, refPoint, figTitle, cLims, ...
                                                   xLimsPlot, yLimsPlot, zLimsPlot, normDims, figSave);
        end
        
    end
    clear i;
    
else
    disp('Skipping Mode Presentation');
end


%% Local Functions

function modes = inputModes(Nt)

    modes = str2num(input('    Input Desired Modes [Row Vector Form]: ', 's')); %#ok<ST2NM>
    
    if isempty(modes) || any(isnan(modes)) || ~isrow(modes) > 1 || any(modes <= 0) || any(modes > Nt)
        disp('        WARNING: Invalid Entry');
        
        modes = -1;
    end

end