%% Planar Velocity POD Calculator v3.0
% ----
% Perform Proper Orthogonal Decomposition on Previously Processed Planar Velocity Fields
% (Collected Using 'planarVelocityProcessing')


%% Preamble

run preamble;

%#ok<*UNRCH>

flipMode = false; % Present Both Orientations of Mode(s)

viewComps = false; % Present Individual Component Views

normDims = true; % Normalise Spatial Dimensions in Plots

figSave = false; % Save .fig File(s)

disp('===================================');
disp('Planar Pressure POD Calculator v3.0');
disp('===================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats
% v3.0 - Offloaded Reconstruction and Corrected Inconsistent Normalisation


%% Acquire Velocity Data

disp('Velocity Data Acquisition');
disp('--------------------------');

valid = false;
while ~valid
    disp(' ');

    [fileName, filePath] = uigetfile([saveLoc, '/Numerical/MATLAB/planarVelocityFields/*.mat'], ...
                                     'Select Velocity Data');

    if contains(filePath, '/Numerical/MATLAB/planarVelocityFields/')
        disp(['Loading ''', fileName, '''...']);

        campaignID = load([filePath, fileName], 'campaignID').campaignID;
        caseID = load([filePath, fileName], 'caseID').caseID;
        dataID = load([filePath, fileName], 'dataID').dataID;
        planeID = load([filePath, fileName], 'planeID').planeID;
        uData = load([filePath, fileName], 'uData').uData;
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

if contains(planeID, 'X_')
    orientation = 'YZ';
elseif contains(planeID, 'Y_')
    orientation = 'XZ';
else
    orientation = 'XY';
end

switch orientation
    
    case 'YZ'
        xLimsData = uData.(planeID).positionGrid(1,1);
        yLimsData = [min(uData.(planeID).positionGrid(:,2)); max(uData.(planeID).positionGrid(:,2))];
        zLimsData = [min(uData.(planeID).positionGrid(:,3)); max(uData.(planeID).positionGrid(:,3))];
        
    case 'XZ'
        xLimsData = [min(uData.(planeID).positionGrid(:,1)); max(uData.(planeID).positionGrid(:,1))];
        yLimsData = uData.(planeID).positionGrid(1,2);
        zLimsData = [min(uData.(planeID).positionGrid(:,3)); max(uData.(planeID).positionGrid(:,3))];
        
    case 'XY'
        xLimsData = [min(uData.(planeID).positionGrid(:,1)); max(uData.(planeID).positionGrid(:,1))];
        yLimsData = [min(uData.(planeID).positionGrid(:,2)); max(uData.(planeID).positionGrid(:,2))];
        zLimsData = uData.(planeID).positionGrid(1,3);
        
end

disp(' ');
disp(' ');


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry(geoLoc);

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

% Initialise POD Variables
PODdata = uData.(planeID); clear uData;
cellSize = cellSize.(planeID);

disp(' ');

% Perform Planar Snapshot POD
[fig, PODdata, modesEnergetic, modes80percent, Ns, Nt] = performPOD(fig, PODdata, {'u', 'v', 'w'}, ...
                                                                    'vector', planeID, figSave);

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
        
        if ~exist([saveLoc, '/Numerical/MATLAB/planarVelocityPOD/', campaignID, '/', caseID, '/', dataID], 'dir')
            mkdir([saveLoc, '/Numerical/MATLAB/planarVelocityPOD/', campaignID, '/', caseID, '/', dataID]);
        end
        
        disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/planarVelocityPOD/', campaignID, '/', caseID, '/', dataID, '/', planeID, '.mat']);
        save([saveLoc, '/Numerical/MATLAB/planarVelocityPOD/', campaignID, '/', caseID, '/', dataID, '/', planeID, '.mat'], ...
              'campaignID', 'caseID', 'dataID', 'planeID', 'PODdata', 'cellSize', 'sampleInt', 'timePrecision', '-v7.3', '-noCompression');
        disp('        Success');
        
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
        
        parts = fieldnames(geometry);
        for i = 1:height(parts)
            geometry.(parts{i}).vertices = geometry.(parts{i}).vertices / normLength;
        end
        clear i parts;

        xDims = xDims / normLength;
        yDims = yDims / normLength;
        zDims = zDims / normLength;

        xLimsData = xLimsData / normLength;
        yLimsData = yLimsData / normLength;
        zLimsData = zLimsData / normLength;

        cellSize.target = cellSize.target / normLength;
        cellSize.y = cellSize.y / normLength;
        cellSize.z = cellSize.z / normLength;
        cellSize.area = cellSize.area / (normLength^2);

        PODdata.positionGrid = PODdata.positionGrid / normLength;
    end
    
end

disp(' ');
disp(' ');


%% Present POD Modes

disp('Mode Presentation');
disp('------------------');

disp(' ');

if plotModes
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
    
    nComponents = 1;
    
    switch orientation
        
        case 'YZ'
            component = 'u';
            
        case 'XZ'
            component = 'v';
            
        case 'XY'
            component = 'w';
            
    end
    
    mapPerim = [];
    nPlanes = 1;
    planeNo = 1;
    cMap = cool2warm(32);
    streamlines = true;
    cLims = [-1; 1];
    
    switch orientation
        
        case 'YZ'
            xLimsPlot = [0.3; 4.625766283524905];
            yLimsPlot = [-0.5; 0.5];
            zLimsPlot = [0; 0.5];
            
        case {'XZ', 'XY'}
            xLimsPlot = [0.3; 1.2];
            yLimsPlot = [-0.3; 0.3];
            zLimsPlot = [0; 0.5];
                
    end
    
    if ~normDims
        xLimsPlot = xLimsPlot * normLength;
        yLimsPlot = yLimsPlot * normLength;
        zLimsPlot = zLimsPlot * normLength;
    end
    
    for i = nModes
        disp(['    Presenting Mode #', num2str(i), '...']);
        
        vectorData = [PODdata.POD.phi((1:Ns),i), ...
                      PODdata.POD.phi(((Ns + 1):(2 * Ns)),i), ...
                      PODdata.POD.phi((((2 * Ns) + 1):end),i)];
        vectorData = vectorData / max(abs(vectorData(:)));
        figName = [planeID, '_POD_uvw_M', num2str(i)];
        figTitle = [num2str(round(PODdata.POD.modeEnergy(i), 2), '%.2f'), '%'];
        
        % Three-Component Plots
        [fig, planeNo] = plotPlanarVectorField(orientation, positionData, vectorData, spatialRes, ...
                                               xLimsData, yLimsData, zLimsData, nComponents, component, ...
                                               mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, ...
                                               streamlines, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                               zLimsPlot, normDims, figSave);
        
        if flipMode
            [fig, planeNo] = plotPlanarVectorField(orientation, positionData, vectorData, spatialRes, ...
                                                   xLimsData, yLimsData, zLimsData, nComponents, component, ...
                                                   mapPerim, nPlanes, planeNo, fig, [figName, '_Flip'], flipud(cMap), geometry, ...
                                                   streamlines, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                                   zLimsPlot, normDims, figSave);
        end
        
        % Single-Component Plots
        if viewComps
            
            for j = ['u', 'v', 'w']
                figName = [planeID, '_POD_', j, '_M', num2str(i)];

                [fig, planeNo] = plotPlanarVectorField(orientation, positionData, vectorData, spatialRes, ...
                                                       xLimsData, yLimsData, zLimsData, nComponents, j, ...
                                                       mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, ...
                                                       false, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                                       zLimsPlot, normDims, figSave);

                if flipMode
                    [fig, planeNo] = plotPlanarVectorField(orientation, positionData, vectorData, spatialRes, ...
                                           xLimsData, yLimsData, zLimsData, nComponents, j, ...
                                           mapPerim, nPlanes, planeNo, fig, [figName, '_Flip'], flipud(cMap), geometry, ...
                                           false, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                           zLimsPlot, normDims, figSave);
                end

            end
            clear j;
            
        end
        
    end
    clear i;
    
    disp(' ');
end
    
if ~plotModes
    disp('Skipping Mode Presentation');
    
    disp(' ');
end


%% Local Functions

function modes = inputModes(Nt)

    modes = str2num(input('    Input Desired Modes [Row Vector Form]: ', 's')); %#ok<ST2NM>
    
    if isempty(modes) || any(isnan(modes)) || ~isrow(modes) > 1 || any(modes <= 0) || any(modes > Nt)
        disp('        WARNING: Invalid Entry');
        
        modes = -1;
    end

end