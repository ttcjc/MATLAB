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


%% Planar Lagrangian Contaminant POD Calculator v2.0

figSave = false; % Save .fig File(s)

disp('======================================');
disp('Planar Contaminant POD Calculator v2.0');
disp('======================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit (Base Contamination and Far-Field Extraction Plane)
% v1.1 - Rename and Restructure to Account for Changes to 'mapData' Output
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats


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

disp('Contaminant Map Acquisition');
disp('----------------------------');

switch format
    
    case {'A', 'B'}
        
        valid = false;
        while ~valid
            disp(' ');
            [fileName, filePath] = uigetfile([saveLocation, '/Numerical/MATLAB/planarContaminantMap/*.mat'], ...
                                             'Select Map Data');

            switch format

                case 'A'

                    if contains(filePath, '/base')
                        disp(['Loading ''', fileName, '''...']);
                        
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
            [fileName, filePath] = uigetfile([saveLocation, '/Experimental/MATLAB/planarContaminantMap/*.mat'], ...
                                             'Select Map Data');
            
            if contains(filePath, 'Hz')
                disp(['Loading ''', fileName, '''...']);
                
                caseID = load([filePath, fileName], 'caseID').caseID;
                planeID = load([filePath, fileName], 'planeID').planeID;
                dataID = load([filePath, fileName], 'dataID').dataID;
                mapData = load([filePath, fileName], 'mapData').mapData;
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

disp('POD Options');
disp('------------');

% Select Variable of Interest
PODvar = fieldnames(mapData.prime);

valid = false;
while ~valid
    disp(' ');
    [index, valid] = listdlg('listSize', [300, 300], ...
                             'selectionMode', 'single', ...
                             'name', 'Select Variable for Decomposition', ...
                             'listString', PODvar);

    if ~valid
        disp('WARNING: No Mapping Variable Selected');
    end
    
end
clear valid;

PODvar = PODvar{index};

disp(['Variable of Interest: ', PODvar]);

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

        xLimsData = mapData.positionGrid(1,1);
        yLimsData = [min(mapData.positionGrid(:,2)); max(mapData.positionGrid(:,2))];
        zLimsData = [min(mapData.positionGrid(:,3)); max(mapData.positionGrid(:,3))];
        
end

% Initialise POD Variables
PODdata.positionGrid = mapData.positionGrid;
PODdata.time = mapData.time;
PODdata.(PODvar).inst = mapData.inst.(PODvar);
PODdata.(PODvar).mean = mapData.mean.(PODvar);
PODdata.(PODvar).prime = mapData.prime.(PODvar);

clear mapData;

disp(' ');

% Perform Planar Snapshot POD
switch format

    case 'A'
        [fig, PODdata, modesEnergetic, modes80percent, Ns, Nt] = performPOD(fig, PODdata, PODvar, ...
                                                                            'scalar', 'Base', figSave);

    case 'B'
        [fig, PODdata, modesEnergetic, modes80percent, Ns, Nt] = performPOD(fig, PODdata, PODvar, ...
                                                                            'scalar', planeID, figSave);
    
    case 'C'
        [fig, PODdata, modesEnergetic, modes80percent, Ns, Nt] = performPOD(fig, PODdata, PODvar, ...
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
    % Define Plot Limits
    switch format

        case 'A'
            orientation = 'YZ';

            if contains(caseID, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 1.43075];
                yLimsPlot = [-0.2445; 0.2445];
                zLimsPlot = [0; 0.389];
            end

        case 'B'
            orientation = 'YZ';

            if contains(caseID, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 2.73575];
                yLimsPlot = [-0.522; 0.522];
                zLimsPlot = [0; 0.6264];
            end
            
        case 'C'
            orientation = 'YZ';
            
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
    figTitle = '.'; % Leave Blank ('-') for Formatting Purposes
    cLims = [-1; 1];

    for i = nModes
        disp(['    Presenting Mode #', num2str(i), '...']);

        scalarData = rescale(PODdata.phi_mode(:,i), -1, 1);

        switch format

            case 'A'
                figName = ['Base_POD_', PODvar, '_M', num2str(i)];

            case 'B'
                figName = [planeID, '_POD_', PODvar, '_M', num2str(i)];
            
            case 'C'
                figName = [planeID, '_POD_', PODvar, '_M', num2str(i)];

        end

        figSubtitle = [num2str(round(PODdata.modeEnergy(i), 2), '%.2f'), '%'];

        [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                               xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                               planeNo, fig, figName, cMap, geometry, contourlines, ...
                                               refPoint, figTitle, figSubtitle, cLims, xLimsPlot, ...
                                               yLimsPlot, zLimsPlot, normDims, figSave);
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
                
                if ~exist([saveLocation, '/Numerical/MATLAB/planarContaminantPOD/', caseID, '/base/', PODvar], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/planarContaminantPOD/', caseID, '/base/', PODvar]);
                end
                
            case 'B'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/planarContaminantPOD/', caseID, '/', planeID, '/', PODvar], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/planarContaminantPOD/', caseID, '/', planeID, '/', PODvar]);
                end
                
            case 'C'
                
                if ~exist([saveLocation, '/Experimental/MATLAB/planarContaminantPOD/', caseID, '/', PODvar], 'dir')
                    mkdir([saveLocation, '/Experimental/MATLAB/planarContaminantPOD/', caseID, '/', PODvar]);
                end
                
        end
        
        switch format
            
            case 'A'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/planarContaminantPOD/', caseID, '/base/', PODvar, '/', dataID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/planarContaminantPOD/', caseID, '/base/', PODvar, '/', dataID, '.mat'], ...
                      'caseID', 'dataID', 'PODdata', 'sampleInt', 'timePrecision', 'normDims', '-v7.3', '-noCompression');
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/planarContaminantPOD/', caseID, '/', planeID, '/', PODvar, '/', dataID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/planarContaminantPOD/', caseID, '/', planeID, '/', PODvar, '/', dataID, '.mat'], ...
                      'caseID', 'planeID', 'dataID', 'PODdata', 'sampleInt', 'timePrecision', 'normDims', '-v7.3', '-noCompression');
                disp('        Success');
            
            case 'C'
                disp(['    Saving to: ', saveLocation, '/Experimental/MATLAB/planarContaminantPOD/', caseID, '/', PODvar, '/', dataID, '.mat']);
                save([saveLocation, '/Experimental/MATLAB/planarContaminantPOD/', caseID, '/', PODvar, '/', dataID, '.mat'], ...
                      'caseID', 'planeID', 'dataID', 'PODdata', 'sampleInt', 'timePrecision', 'normDims', '-v7.3', '-noCompression');
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


%% Select Reconstruction Options

disp('Reconstruction Options');
disp('-----------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Perform Field Reconstruction Using N Modes? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        return;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
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


%% Perform Field Reconstruction

disp('Field Reconstruction');
disp('---------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

disp(' ');

disp('    Initialising...');

% Initialise Reconstruction Variables
reconData.positionGrid = PODdata.positionGrid;
reconData.time = PODdata.time;
reconData.(PODvar).inst = cell(Nt,1);

for i = 1:Nt
    reconData.(PODvar).inst{i} = PODdata.(PODvar).mean;
end
clear i;

disp(' ');

% Perform Field Reconstruction
if width(nModes) <= 5
    reconData = reconstructPOD(reconData, PODdata, PODvar, nModes, Ns, Nt, 'scalar', true);
else
    reconData = reconstructPOD(reconData, PODdata, PODvar, nModes, Ns, Nt, 'scalar', false);
end

disp(' ');

% Calculate Calculate Reconstructed Time Average
disp('    Calculating Reconstructed Time Average');

reconData.(PODvar).mean = zeros([Ns,1]);

for i = 1:Nt
    reconData.(PODvar).mean = reconData.(PODvar).mean + reconData.(PODvar).inst{i};
end
clear i;

reconData.(PODvar).mean = reconData.(PODvar).mean / Nt;

disp(' ');

% Calculate Instantaneous Spray Density Fluctuations
disp('    Calculating Reconstructed Instantaneous Fluctuations');

reconData.(PODvar).prime = reconData.(PODvar).inst;

for i = 1:Nt
    reconData.(PODvar).prime{i} = reconData.(PODvar).inst{i} - reconData.(PODvar).mean;
end
clear i;

disp(' ');

% Calculate RMS of Reconstructed Field
disp('    Calculating RMS of Reconstructed Field');

reconData.(PODvar).RMS = zeros([Ns,1]);

for i = 1:Nt
    reconData.(PODvar).RMS  = reconData.(PODvar).RMS + reconData.(PODvar).prime{i}.^2;
end
clear i;

reconData.(PODvar).RMS = sqrt((1 / Nt) * reconData.(PODvar).RMS);

if (any(strcmp(format, {'A', 'B'})) && any(strcmp(PODvar, {'mass', 'massNorm'}))) || strcmp(format, 'C')
    disp(' ');

    disp('    Calculating Reconstructed Centre of Mass...');

    % Initialise Progress Bar
    wB = waitbar(0, 'Calculating Reconstructed Centre of Mass', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    
    % Calculate Reconstructed CoM
    reconData.CoM.inst = cell(Nt,1); reconData.CoM.inst(:) = {zeros([1,3])};
    
    for i = 1:Nt
        reconData.CoM.inst{i}(1) = reconData.positionGrid(1,1);
        reconData.CoM.inst{i}(2) = sum(reconData.(PODvar).inst{i} .* reconData.positionGrid(:,2)) / sum(reconData.(PODvar).inst{i});
        reconData.CoM.inst{i}(3) = sum(reconData.(PODvar).inst{i} .* reconData.positionGrid(:,3)) / sum(reconData.(PODvar).inst{i});
        
        % Update Waitbar
        waitbar((i / Nt), wB);
    end
    clear i;
    
    delete(wB);
    
    disp(' ');

    % Calculate Time-Averaged Centre of Spray
    disp('    Calculating Time-Averaged Centre of Spray');

    reconData.CoM.mean = zeros([1,3]);

    reconData.CoM.mean(1) = mapData.positionGrid(1,1);
    reconData.CoM.mean(2) = sum(reconData.(PODvar).mean .* reconData.positionGrid(:,2)) / sum(reconData.(PODvar).mean);
    reconData.CoM.mean(3) = sum(reconData.(PODvar).mean .* reconData.positionGrid(:,3)) / sum(reconData.(PODvar).mean);
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
    % Define Plot Limits
    switch format

        case 'A'
            orientation = 'YZ';

            if contains(caseID, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 1.43075];
                yLimsPlot = [-0.2445; 0.2445];
                zLimsPlot = [0; 0.389];
            end

        case 'B'
            orientation = 'YZ';

            if contains(caseID, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 2.73575];
                yLimsPlot = [-0.522; 0.522];
                zLimsPlot = [0; 0.6264];
            end
            
        case 'C'
            orientation = 'YZ';
            
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
    cMap = flipud(viridis(32));
    refPoint = [];
    figTitle = '.'; % Leave Blank ('-') for Formatting Purposes
    
    % Remove Geometry From Empty Tunnel
    if contains(caseID, 'ET')
        geometry = [];
    end

end

if plotMean
    disp('    Presenting Average Reconstructed  Field...');

    scalarData = reconData.(PODvar).mean;
    figName = ['Recon_Average_', caseID];
    figSubtitle = ' ';
    
    switch format
        
        case {'A', 'B'}
            cLims = [0; max(cellfun(@max, reconData.(PODvar).inst))];
            contourlines = [];
            
        case 'C'
            cLims = [0; 1.05];
            contourlines = [0.02; 0.02];
            
    end

    [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                           xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                           planeNo, fig, figName, cMap, geometry, contourlines, ...
                                           refPoint, figTitle, figSubtitle, cLims, xLimsPlot, ...
                                           yLimsPlot, zLimsPlot, normDims, figSave);
    
    disp(' ');
end

if plotRMS
    disp('    Presenting RMS of Reconstructed Field...');
    
    scalarData = reconData.(PODvar).RMS;
    figName = ['Recon_RMS_', caseID];
    contourlines = [];
    figSubtitle = ' ';
    cLims = 'auto';

    [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                           xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                           planeNo, fig, figName, cMap, geometry, contourlines, ...
                                           refPoint, figTitle, figSubtitle, cLims, xLimsPlot, ...
                                           yLimsPlot, zLimsPlot, normDims, figSave);
    
    disp(' ');
end

if plotInst
    disp('    Presenting Instantaneous Reconstructed Field...');
    
    switch format
        
        case {'A', 'B'}
            cLims = [0; max(cellfun(@max, reconData.(PODvar).inst))];
            contourlines = [];
            
        case 'C'
            cLims = [0; 3.6];
            contourlines = [0.02; 0.02];
            
    end

    figHold = fig;

    for i = startFrame:endFrame
        
        if i ~= startFrame
            clf(fig);
            fig = figHold;
        end
        
        scalarData = reconData.(PODvar).inst{i};
        figTime = num2str(reconData.time(i), ['%.', num2str(timePrecision), 'f']);
        
        switch format
            
            case 'A'
                figName = ['Base_', PODvar, '_Reconstruction_T', erase(figTime, '.')];
            
            case 'B'
                figName = [planeID, '_', PODvar, '_Reconstruction_T', erase(figTime, '.')];
                
            case 'C'
                figName = [planeID, '_', PODvar, '_Reconstruction_T', erase(figTime, '.')];
        
        end
        
        figSubtitle = [figTime, ' \it{s}'];
        
        [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                               xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                               planeNo, fig, figName, cMap, geometry, contourlines, ...
                                               refPoint, figTitle, figSubtitle, cLims, xLimsPlot, ...
                                               yLimsPlot, zLimsPlot, normDims, figSave);
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
                
                if ~exist([saveLocation, '/Numerical/MATLAB/planarContaminantReconstruction/', caseID, '/base/', PODvar], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/planarContaminantReconstruction/', caseID, '/base/', PODvar]);
                end
                
            case 'B'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/planarContaminantReconstruction/', caseID, '/', planeID, '/', PODvar], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/planarContaminantReconstruction/', caseID, '/', planeID, '/', PODvar]);
                end
                
            case 'C'
                
                if ~exist([saveLocation, '/Experimental/MATLAB/planarContaminantReconstruction/', caseID, '/', PODvar], 'dir')
                    mkdir([saveLocation, '/Experimental/MATLAB/planarContaminantReconstruction/', caseID, '/', PODvar]);
                end
                
        end
        
        switch format
            
            case 'A'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/planarContaminantReconstruction/', caseID, '/base/', PODvar, '/', dataID, '_', mat2str(nModes), '.mat']);
                save([saveLocation, '/Numerical/MATLAB/planarContaminantReconstruction/', caseID, '/base/', PODvar, '/', dataID, '_', mat2str(nModes), '.mat'], ...
                     'caseID', 'dataID', 'reconData', 'nModes', 'sampleInt', 'dLims', 'timePrecision', 'normDims', '-v7.3', '-noCompression');
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/planarContaminantReconstruction/', caseID, '/', planeID, '/', PODvar, '/', dataID, '_', mat2str(nModes), '.mat']);
                save([saveLocation, '/Numerical/MATLAB/planarContaminantReconstruction/', caseID, '/', planeID, '/', PODvar, '/', dataID, '_', mat2str(nModes), '.mat'], ...
                     'caseID', 'planeID', 'dataID', 'reconData', 'nModes', 'sampleInt', 'dLims', 'timePrecision', 'normDims', '-v7.3', '-noCompression');
                disp('        Success');
                
            case 'C'
                disp(['    Saving to: ', saveLocation, '/Experimental/MATLAB/planarContaminantReconstruction/', caseID, '/', PODvar, '/', dataID, '_', mat2str(nModes), '.mat']);
                save([saveLocation, '/Experimental/MATLAB/planarContaminantPOD/', caseID, '/', PODvar, '/', dataID, '_', mat2str(nModes), '.mat'], ...
                     'caseID', 'planeID', 'dataID', 'reconData', 'nModes', 'sampleInt', 'timePrecision', 'normDims', '-v7.3', '-noCompression');
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