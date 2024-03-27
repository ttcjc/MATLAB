%% Lagrangian Line of Sight Obstruction Calculator v2.3
% ----
% Calculate Mass of Spray Contained Along Rays Projecting From a Point of Interest to an Upstream Plane
% (Based on Volumetric Spray Data Collected Using 'volumeFieldGenerator')


%% Preamble

run preamble;

samplesPerCell = 16; % Used During Numerical Integration

normDims = true; % Normalise Spatial Dimensions in Plots

figSave = false; % Save .fig File(s)

disp('=========================================');
disp('Line of Sight Obstruction Calculator v2.2');
disp('=========================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Added Support for Full-Scale Windsor Model Simulations
% v2.1 - Minor Update to Shift Preamble Into Separate Script
% v2.2 - Update To Correct Inconsistent Normalisation Throughout Repository
% V2.3 - Offloaded Rays of Interest Calculations


%% Select Region of Interest

disp('Region of Interest');
disp('-------------------');

disp(' ');

disp('Possible Regions of Interest:');
disp('    A: Near Wake');
disp('    B: Mid Wake');
disp('    C: Far Wake (Full-Scale Only)');

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Select Region of Interest [A/B]: ', 's');

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


%% Acquire Volume Field

disp('Volume Field Acquisition');
disp('-------------------------');

valid = false;
while ~valid
    disp(' ');
    [fileName, filePath] = uigetfile([saveLoc, '/Numerical/MATLAB/volumeField/*.mat'], ...
                                      'Select Volumetric Data');
    
    switch format
        
        case 'A'
            
            if contains(filePath, '/nearWake')
                valid = true;
            else
                disp('WARNING: Invalid File Selection');
                clear fileName filePath;
            end
            
        case 'B'
            
            if contains(filePath, '/midWake')
                valid = true;
            else
                disp('WARNING: Invalid File Selection');
                clear fileName filePath;
            end
            
        case 'C'
            
            if contains(filePath, '/farWake')
                valid = true;
            else
                disp('WARNING: Invalid File Selection');
                clear fileName filePath;
            end
            
    end
    
    disp(['Loading ''', fileName, '''...']);
    
    campaignID = load([filePath, fileName], 'campaignID').campaignID;
    caseID = load([filePath, fileName], 'caseID').caseID;
    dataID = load([filePath, fileName], 'dataID').dataID;
    volumeData = load([filePath, fileName], 'volumeData').volumeData;
    cellSize = load([filePath, fileName], 'cellSize').cellSize;
    sampleInt = load([filePath, fileName], 'sampleInt').sampleInt;
    timePrecision = load([filePath, fileName], 'timePrecision').timePrecision;
    dLims = load([filePath, fileName], 'dLims').dLims;
    
    disp('    Success');
end
clear valid;

disp(' ');
disp(' ');


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry(geoLoc);

disp(' ');
disp(' ');


%% Select Plane of Interest

% Select Plane of Interest
obstructData.PoV.targetPlane = identifyVolumeSlices(volumeData.positionGrid, spacePrecision, false);
planeID = fieldnames(obstructData.PoV.targetPlane); planeID = planeID{1};
obstructData.PoV.targetPlane = obstructData.PoV.targetPlane.(planeID);
clear planeID;

% Extract Planar Position Data
orientation = obstructData.PoV.targetPlane.orientation;

switch orientation
    
    case 'YZ'
        index = find(volumeData.positionGrid(:,1) == obstructData.PoV.targetPlane.position);
        
    case 'XZ'
        index = find(volumeData.positionGrid(:,2) == obstructData.PoV.targetPlane.position);
        
    case 'XY'
        index = find(volumeData.positionGrid(:,3) == obstructData.PoV.targetPlane.position);
        
end

obstructData.positionGrid = volumeData.positionGrid(index,:);

nCells = height(obstructData.positionGrid);

disp(' ');
disp(' ');


%% Select Origin Point

disp('Ray Origin Definition');
disp('----------------------');

% Select Ray Origin Point
obstructData.PoV.originPoint = zeros([1,3]);

valid = false;
while ~valid
    disp(' ')
    disp('Specify Ray Origin Point:')
    
    obstructData.PoV.originPoint(1) = inputPos('X');
    obstructData.PoV.originPoint(2) = inputPos('Y');
    obstructData.PoV.originPoint(3) = inputPos('Z');

    if (obstructData.PoV.originPoint(1) < min(volumeData.positionGrid(:,1)) || ...
        obstructData.PoV.originPoint(1) > max(volumeData.positionGrid(:,1))) || ...
       (obstructData.PoV.originPoint(2) < min(volumeData.positionGrid(:,2)) || ...
        obstructData.PoV.originPoint(2) > max(volumeData.positionGrid(:,1))) || ...
       (obstructData.PoV.originPoint(3) < min(volumeData.positionGrid(:,3)) || ...
        obstructData.PoV.originPoint(3) > max(volumeData.positionGrid(:,1)))
        disp('        WARNING: Origin Point Lies Outside Volume');
        
        continue;
    end
    
    switch orientation
        
        case 'YZ'
            
            if obstructData.PoV.originPoint(1) < obstructData.PoV.targetPlane.position
                disp('        WARNING: Origin Point Must Lie in the Positive X-Direction of the Target Plane');
                
                continue;
            end
            
        case 'XZ'
            
            if obstructData.PoV.originPoint(2) > obstructData.PoV.targetPlane.position
                disp('        WARNING: Origin Point Must Lie in the Negative Y-Direction of the Target Plane');
                
                continue;
            end
            
        case 'XY'
            
            if obstructData.PoV.originPoint(3) > obstructData.PoV.targetPlane.position
                disp('        WARNING: Origin Point Must Lie in the Positive Z-Direction of the Target Plane');
                
                continue;
            end
            
    end
    
    valid = true;    
end
clear valid;

disp(' ');
disp(' ');


%% Calculate Line of Sight Obstruction

disp('Line of Slight Obstruction Calculation');
disp('---------------------------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

maxNumCompThreads(nProc);

evalc('parpool(''threads'');');

%%%%

disp(' ');

disp('    Initialising...');

nTimes = height(volumeData.time);

dL = cellSize.volume^(1 / 3) / samplesPerCell;

% Check if Plane of Interest Intersects Geometry
removeIntersect = false;

switch orientation
    
    case 'YZ'
        
        if obstructData.PoV.targetPlane.position <= xDims(2)
            removeIntersect = true;
        end
        
    case 'XZ'
        
        if obstructData.PoV.targetPlane.position >= yDims(1)
            removeIntersect = true;
        end
        
    case 'XY'
        
        if obstructData.PoV.targetPlane.position <= zDims(2)
            removeIntersect = true;
        end
        
end

% Remove Erroneous Data From Cells Intersecting Geometry
if removeIntersect
    disp('        Removing Erroneous Data From Grid Cells Intersecting Geometry...');

    % Perform Removal
    volumeDataVars = fieldnames(volumeData);
    nonFieldVars = {'positionGrid'; 'time'};
    fieldVars = setdiff(volumeDataVars, nonFieldVars);
    clear volumeDataVars nonFieldVars;
    
    parts = fieldnames(geometry);
    for i = 1:height(parts)
        DT = delaunay(geometry.(parts{i}).vertices);

        index = ~isnan(tsearchn(geoPoints, DT, volumeData.positionGrid));

        for j = 1:height(fields)
            volumeData.(fieldVars{j}).mean(index,:) = NaN;

            for k = 1:nTimes
                volumeData.(fieldVars{j}).inst{k}(index,:) = NaN;
            end

        end
        clear j;

    end
    clear i parts;
    
    clear fieldVars;    
end

% Update Map Boundaries
switch orientation

    case 'YZ'
        xLimsData = obstructData.PoV.targetPlane.position;
        yLimsData = [min(obstructData.positionGrid(:,2)); max(obstructData.positionGrid(:,2))];
        zLimsData = [min(obstructData.positionGrid(:,3)); max(obstructData.positionGrid(:,3))];
        
    case 'XZ'
        xLimsData = [min(obstructData.positionGrid(:,1)); max(obstructData.positionGrid(:,1))];
        yLimsData = obstructData.PoV.targetPlane.position;
        zLimsData = [min(obstructData.positionGrid(:,3)); max(obstructData.positionGrid(:,3))];
        
    case 'XY'
        xLimsData = [min(obstructData.positionGrid(:,1)); max(obstructData.positionGrid(:,1))];
        yLimsData = [min(obstructData.positionGrid(:,2)); max(obstructData.positionGrid(:,2))];
        zLimsData = obstructData.PoV.targetPlane.position;
        
end

disp (' ');

% Reshape Position Data for Improved Interpolation Performance
disp('    Reshaping Position Data for Improved Interpolation Performance...');

gridShape = [height(unique(volumeData.positionGrid(:,1))), ...
             height(unique(volumeData.positionGrid(:,2))), ...
             height(unique(volumeData.positionGrid(:,3)))];

x = reshape(volumeData.positionGrid(:,1), gridShape);
y = reshape(volumeData.positionGrid(:,2), gridShape);
z = reshape(volumeData.positionGrid(:,3), gridShape);

disp(' ');

% Calculate Instantaneous Line of Sight
disp('    Calculating Instantaneous Obstruction...');

obstructData.time = volumeData.time;

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Obstruction', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);
parforWaitBar(wB, nTimes);

% Perform Calculation
density = cell(nTimes,1); density(:) = {zeros([nCells,1])};

densityVolume = volumeData.density.inst;
positionGrid = obstructData.positionGrid;
pointPosition = obstructData.PoV.originPoint;
parfor i = 1:nTimes
    densityField = reshape(full(densityVolume{i}), gridShape);
    densityInterp = griddedInterpolant(x, y, z, densityField, 'linear', 'none');
    
    for j = 1:height(positionGrid)
        dirVec = positionGrid(j,:) - pointPosition;
        distFull = sqrt(dirVec(1)^2 + dirVec(2)^2 + dirVec(3)^2);

        dist = (dL:dL:distFull)';
        samplePoints = pointPosition + (dist * (dirVec / distFull));
        
        density{i}(j) = sum((densityInterp(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3))) * dL);
    end
    
    % Remove Unnecessary Data
    densityVolume{i} = [];
    
    % Update Waitbar
    send(dQ, []);
end
clear densityVolume positionGrid pointPosition;

delete(wB);

clear volumeData;

obstructData.density.inst = density; clear density;

disp(' ');

disp('    Calculating Time-Averaged Obstruction...');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Time-Averaged Obstruction', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Perform Calculation
obstructData.density.mean = zeros([nCells,1]);

for i = 1:nTimes
    obstructData.density.mean = obstructData.density.mean + obstructData.density.inst{i};
    
    % Update Waitbar
    waitbar((i / nTimes), wB);
end
clear i;

delete(wB);

obstructData.density.mean = obstructData.density.mean / nTimes;

disp(' ');

% Generate Time-Averaged Spray Map
disp('    Generating Fluctuating Obstruction...');

% Calculate Instantaneous Field Fluctuations
disp('        Calculating Instantaneous Obstruction Fluctuations');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Obstruction Fluctuations', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Perform Calculation
obstructData.density.prime = obstructData.density.inst;

for i = 1:nTimes
    obstructData.density.prime{i} = obstructData.density.prime{i} - obstructData.density.mean;
    
    % Update Waitbar
    waitbar((i / nTimes), wB);
end
clear i;

delete(wB);

% Calculate RMS of Field Variables
disp('        Calculating RMS of Obstruction Fluctuations');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating RMS of Obstruction Fluctuations', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Perform Calculation
obstructData.density.RMS = zeros([nCells,1]);

for i = 1:nTimes
    obstructData.density.RMS = obstructData.density.RMS + obstructData.density.prime{i}.^2;
    
    % Update Waitbar
    waitbar((i / nTimes), wB);
end
clear i;

delete(wB);

obstructData.density.RMS = sqrt((1 / nTimes) * obstructData.density.RMS);

switch orientation
    
    case 'YZ'
        obstructID = ['target_X' ...
                      replace(num2str(round(obstructData.PoV.targetPlane.position, 3)), '.', '_'), ...
                      '_origin_X', replace(num2str(round(obstructData.PoV.originPoint(1), 3)), '.', '_'), ...
                      '_Y', replace(num2str(round(obstructData.PoV.originPoint(2), 3)), '.', '_'), ...
                      '_Z', replace(num2str(round(obstructData.PoV.originPoint(3), 3)), '.', '_')];
        
    case 'XZ'
        obstructID = ['target_Y' ...
                      replace(num2str(round(obstructData.PoV.targetPlane.position, 3)), '.', '_'), ...
                      '_origin_X', replace(num2str(round(obstructData.PoV.originPoint(1), 3)), '.', '_'), ...
                      '_Y', replace(num2str(round(obstructData.PoV.originPoint(2), 3)), '.', '_'), ...
                      '_Z', replace(num2str(round(obstructData.PoV.originPoint(3), 3)), '.', '_')];
        
    case 'XY'
        obstructID = ['target_Z' ...
                      replace(num2str(round(obstructData.PoV.targetPlane.position, 3)), '.', '_'), ...
                      '_origin_X', replace(num2str(round(obstructData.PoV.originPoint(1), 3)), '.', '_'), ...
                      '_Y', replace(num2str(round(obstructData.PoV.originPoint(2), 3)), '.', '_'), ...
                      '_Z', replace(num2str(round(obstructData.PoV.originPoint(3), 3)), '.', '_')];
        
end


obstructData = orderfields(obstructData, {'positionGrid', 'time', 'density', ...
                                          'PoV'});

%%%%

evalc('delete(gcp(''nocreate''));');

executionTime = toc;

disp(' ');

disp(['    Run Time: ', num2str(executionTime), 's']);

disp(' ');

disp('  SUCCESS  ');
disp('***********');

disp(' ');
disp(' ');


%% Save Map Data

disp('Data Save Options');
disp('------------------');

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Save Data for Future Use? [y/n]: ', 's');
    
    if selection == 'n' | selection == 'N' %#ok<OR2>
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        
        % Save Data
        switch format
            
            case 'A'
                
                if ~exist([saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/nearWake/', dataID], 'dir')
                    mkdir([saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/nearWake/', dataID]);
                end
                
                disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/nearWake/', dataID, '/', obstructID, '.mat']);
                
                save([saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/nearWake/', dataID, '/', obstructID, '.mat'], ...
                     'campaignID', 'caseID', 'dataID', 'obstructID', 'obstructData', 'cellSize', 'sampleInt', 'timePrecision', 'dLims', 'dL', '-v7.3', '-noCompression');
                
                disp('        Success');
                
            case 'B'
                
                if ~exist([saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/midWake/', dataID], 'dir')
                    mkdir([saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/midWake/', dataID]);
                end
                
                disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/midWake/', dataID, '/', obstructID, '.mat']);
                
                save([saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/midWake/', dataID, '/', obstructID, '.mat'], ...
                     'campaignID', 'caseID', 'dataID', 'obstructID', 'obstructData', 'cellSize', 'sampleInt', 'timePrecision', 'dLims', 'dL', '-v7.3', '-noCompression');
                
                disp('        Success');
                
            case 'C'
                
                if ~exist([saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/farWake/', dataID], 'dir')
                    mkdir([saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/farWake/', dataID]);
                end
                
                disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/farWake/', dataID, '/', obstructID, '.mat']);
                
                save([saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/farWake/', dataID, '/', obstructID, '.mat'], ...
                     'campaignID', 'caseID', 'dataID', 'obstructID', 'obstructData', 'cellSize', 'sampleInt', 'timePrecision', 'dLims', 'dL', '-v7.3', '-noCompression');
                
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


%% Select Presentation Options

disp('Presentation Options');
disp('---------------------');

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Plot Time-Averaged Obstruction? [y/n]: ', 's');

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
    
    selection = input('Plot RMS of Obstruction? [y/n]: ', 's');

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
    
    selection = input('Plot Instantaneous Obstruction? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotInst = false;
        
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotInst = true;
        
        startFrame = inputFrames(nTimes, 'Start');
        
        if startFrame == -1
            continue;
        end
        
        endFrame = inputFrames(nTimes, 'End');
        
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
        cellSize.volume = cellSize.volume / (normLength^3);
        
        dL = dL / normLength;

        obstructData.positionGrid = obstructData.positionGrid / normLength;
    end
    
end

disp(' ');
disp(' ');


%% Present Obstruction Maps

disp('Map Presentation');
disp('-----------------');

disp(' ');

if plotMean || plotRMS || plotInst    
    orientation = obstructData.PoV.targetPlane.orientation;
    positionData = obstructData.positionGrid;
    
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
    
    mapPerim = [];
    nPlanes = 1;
    planeNo = 1;
    cMap = flipud(viridis(32));
    refPoint = [];
    xLimsPlot = [0.3; 4.6257662];
    yLimsPlot = [-0.5; 0.5];
    zLimsPlot = [0; 0.5];
    
    if ~normDims
        xLimsPlot = xLimsPlot * normLength;
        yLimsPlot = yLimsPlot * normLength;
        zLimsPlot = zLimsPlot * normLength;
    end
    
end

if strcmp(caseID, 'Windsor_SB_wW_Upstream_SC')
    refValue = 1;
elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_uncoupled')
    refValue = 5.064762429899226; % 4L 99%
%     refValue = 14.312439739727616; % 4L 100%
elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_coupled')
    refValue = 8.069761393719082; % 4L 99%
%     refValue = 50.803795281819460; % 4L 100%
elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_halfTread')
    refValue = 4.266535645576186; % 4L 99%
%     refValue = 22.098401361629170; % 4L 100%
elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_20deg')
    refValue = 9.030407587347277; % 4L 99%
%     refValue = 19.479226269422180; % 4L 100%
else
    refValue = prctile(obstructData.density.mean(obstructData.density.mean > 0), 99);
%     refValue = prctile(obstructData.density.mean(obstructData.density.mean > 0), 100);
end

if plotMean
    clear cLims;
    
    disp('    Presenting Time-Averaged Obstruction Map...');
    
    scalarData = obstructData.density.mean;
    figName = ['Average_Obstruction_', obstructID];
    contourlines = [0.02; 0.02] * refValue;
    figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
    
    if strcmp(caseID, 'Windsor_SB_wW_Upstream_SC')
        cLims = 'auto';
    elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_uncoupled')
        cLims = [0; 5.1]; % 4L 99%
%         cLims = [0; 14.4]; % 4L 100%
    elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_coupled')
        cLims = [0; 8.1]; % 4L 99%
%         cLims = [0; 50.8]; % 4L 100%
    elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_halfTread')
        cLims = [0; 4.3]; % 4L 99%
%         cLims = [0; 22.1]; % 4L 100%
    elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_20deg')
        cLims = [0; 9.05]; % 4L 99%
%         cLims = [0; 19.5]; % 4L 100%
    else
        cLims = 'auto';
    end
    
    if ~exist('cLims', 'var')
        cLims = [0; max(scalarData)];
    end
    
    [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                           xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                           planeNo, fig, figName, cMap, geometry, contourlines, ...
                                           refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                           zLimsPlot, normDims, figSave);
    
    disp(' ');
end

if plotRMS
    clear cLims;

    disp('    Presenting RMS of Obstruction Map...');

    scalarData = obstructData.density.RMS / mean(obstructData.density.mean(obstructData.density.mean > 0));
    figName = ['RMS_Obstruction_', obstructID];
    contourlines = [];
    figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes

    if strcmp(caseID, 'Windsor_SB_wW_Upstream_SC')
        cLims = 'auto';
    elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_uncoupled')
        cLims = [0; 5.4];
    elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_coupled')
        cLims = [0; 5.75];
    elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_halfTread')
        cLims = [0; 1.65];
    elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_20deg')
        cLims = [0; 1.4];
    else
        cLims = 'auto';
    end

    if ~exist('cLims', 'var')
        cLims = [0; max(scalarData)];
    end

    [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                           xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                           planeNo, fig, figName, cMap, geometry, contourlines, ...
                                           refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                           zLimsPlot, normDims, figSave);
    
    disp(' ');
end

if plotInst
    clear cLims;

    disp('    Presenting Instantanoues Obstruction Map(s)...');

    contourlines = [];
    
    if strcmp(caseID, 'Windsor_SB_wW_Upstream_SC')
        cLims = 'auto';
    elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_uncoupled')
        cLims = [0; 6.2]; % 4L 99%
%         cLims = [0; 20.9]; % 4L 100%
    elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_coupled')
        cLims = [0; 9.6]; % 4L 99%
%         cLims = [0; 58]; % 4L 100%
    elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_halfTread')
        cLims = [0; 5]; % 4L 99%
%         cLims = [0; 25.25]; % 4L 100%
    elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_20deg')
        cLims = [0; 10.95]; % 4L 99%
%         cLims = [0; 21.9]; % 4L 100%
    else
        cLims = 'auto';
    end

    if ~exist('cLims', 'var')
        instMax = zeros([nTimes,1]);

        for j = 1:nTimes
            instMax(j) = prctile(obstructData.density.inst{j}, 99);
%             instMax(j) = max(obstructData.density.inst{j});
        end
        clear j;
        
        cLims = [0; max(instMax)];
%         cLims = [0; prctile(instMax, 99)];
    end

    figHold = fig;

    for j = startFrame:endFrame

        if j ~= startFrame
            clf(fig);
            fig = figHold;
        end

        scalarData = obstructData.density.inst{j};
        figTime = num2str(obstructData.time(j), ['%.', num2str(timePrecision), 'f']);
        figName = ['Inst_Obstruction_', obstructID, '_T', replace(figTime, '.', '_')];
        figTitle = ['{', figTime, ' \it{s}}'];

        [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                               xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                               planeNo, fig, figName, cMap, geometry, contourlines, ...
                                               refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                               zLimsPlot, normDims, figSave);
    end
    clear j;
    
disp(' ');
end

if ~plotMean && ~plotRMS && ~plotInst
    disp('Skipping Map Presentation...');
    
    disp(' ');
end


%% Local Functions

function pos = inputPos(orientation)

    pos = str2double(input(['    ', orientation, '-Position [m]: '], 's'));
    
    if isnan(pos) || length(pos) > 1
        disp('        WARNING: Invalid Entry');
        
        pos = -1;
    end
    
end


function frameNo = inputFrames(Nt, type)

    frameNo = str2double(input(['    Input Desired ', type, ' Frame [1-', num2str(Nt), ']: '], 's'));
    
    if isnan(frameNo) || frameNo < 1 || frameNo > Nt
        disp('        WARNING: Invalid Entry');
        
        frameNo = -1;
    end

end
