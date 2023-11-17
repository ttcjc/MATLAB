%% Preamble

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

if exist('/mnt/Processing/Data', 'dir')
    saveLoc = '/mnt/Processing/Data';
else
    saveLoc = '~/Data';
end

nProc = maxNumCompThreads - 2; % Number of Processors Used for Process-Based Parallelisation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures


%% Planar Spray POD Calculator v3.0

figSave = false; % Save .fig File(s)

flipMode = true; % Present Both Orientations of Mode(s)

disp('================================');
disp('Planar Spray POD Calculator v3.0');
disp('================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit (Base Contamination and Far-Field Extraction Plane)
% v1.1 - Rename and Restructure to Account for Changes to 'mapData' Output
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats
% v2.1 - Moved Calculation of Instantaneous Variables To Pre-Processing Scripts
% v3.0 - Offloaded Reconstruction To Improve Efficiency


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
                        
                        campaignID = load([filePath, fileName], 'campaignID').campaignID;
                        caseID = load([filePath, fileName], 'caseID').caseID;
                        planeID = load([filePath, fileName], 'planeID').planeID;
                        dataID = load([filePath, fileName], 'dataID').dataID;
                        mapData = load([filePath, fileName], 'mapData').mapData;
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

disp(' ');
disp(' ');


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normDims, normLength] = selectGeometry(normDims);

disp(' ');
disp(' ');


%% Select POD Options

clc

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
        
        mapPerim = boundary(basePoints(:,2), basePoints(:,3), 1);
        mapPerim = basePoints(mapPerim,:);
        basePoly = polyshape(mapPerim(:,2), mapPerim(:,3), 'keepCollinearPoints', true);
        basePoly = polybuffer(basePoly, -0.005, 'jointType', 'square');
        mapPerim = ones([height(basePoly.Vertices),3]) * mapPerim(1,1);
        mapPerim(:,[2,3]) = basePoly.Vertices(:,[1,2]);

        if ~all(mapPerim(1,:) == mapPerim(end,:))
            mapPerim = [mapPerim; mapPerim(1,:)]; % Close Boundary
        end
        
        xLimsData = xDims(2) + 1e-3; % Offset Particles From Base for Better Visibility
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

disp(' ');
disp(' ');


%% Present POD Modes

disp('Mode Presentation');
disp('------------------');

disp(' ');

if plotModes
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

    positionData = PODdata.positionGrid;
    nPlanes = 1;
    planeNo = 1;
    cMap = cool2warm(32);
    contourlines = [];
    refPoint = [];
    cLims = [-1; 1];

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
    disp('    Skipping Mode Presentation');
end

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
                      'campaignID', 'caseID', 'dataID', 'PODdata', 'cellSize', 'sampleInt', 'timePrecision', 'normDims', '-v7.3', '-noCompression');
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/planarSprayPOD/', campaignID, '/', caseID, '/', planeID, '/', field, '/', dataID, '.mat']);
                save([saveLoc, '/Numerical/MATLAB/planarSprayPOD/', campaignID, '/', caseID, '/', planeID, '/', field, '/', dataID, '.mat'], ...
                      'campaignID', 'caseID', 'planeID', 'dataID', 'PODdata', 'cellSize', 'sampleInt', 'timePrecision', 'normDims', '-v7.3', '-noCompression');
                disp('        Success');
            
            case 'C'
                disp(['    Saving to: ', saveLoc, '/Experimental/MATLAB/planarSprayPOD/', campaignID, '/', caseID, '/', field, '/', dataID, '.mat']);
                save([saveLoc, '/Experimental/MATLAB/planarSprayPOD/', campaignID, '/', caseID, '/', field, '/', dataID, '.mat'], ...
                      'campaignID', 'caseID', 'planeID', 'dataID', 'PODdata', 'cellSize', 'sampleInt', 'timePrecision', 'normDims', '-v7.3', '-noCompression');
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

    modes = str2num(input('    Input Desired Modes [Row Vector Form]: ', 's')); %#ok<ST2NM>
    
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
