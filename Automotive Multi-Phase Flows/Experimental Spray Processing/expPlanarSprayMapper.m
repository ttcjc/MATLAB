%% Planar Experimental Spray Mapper v2.1
% ----
% Load, Process and Present Planar Experimental Spray Data


%% Preamble

run preamble;

%#ok<*UNRCH>

normDims = true; % Normalise Spatial Dimensions

normDensity = true; % Normalise Spray Density
    normValue = 0.0052166; % (SB_1.0L_120s_15Hz_02)

figSave = false; % Save .fig File(s);

disp('=====================================');
disp('Planar Experimental Spray Mapper v2.1');
disp('=====================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Update to Support Changes to 'plotPlanarScalarField'
% v2.0 - Update To Support Changes to DaVis Data Format
% v2.1 - Update To Correct Inconsistent Normalisation Throughout Repository


%% Initialise Case

% Select Relevant Geometry and Define Bounding Box
[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry;

disp(' ');
disp(' ');


%% Initialise Experimental Spray Data

if exist('/media/ttcjc/B418689', 'dir')
    dataLoc = '/media/ttcjc/B418689/Data';
else
    dataLoc = '/mnt/Processing/Data';
end

maxNumCompThreads(nProc);

[campaignID, caseID, planeID, ...
 expSprayData, sampleFreq] = initialiseExpSprayData(saveLoc, dataLoc, maxNumCompThreads);

disp(' ');
disp(' ');


%% Select Mapping Options

disp('Mapping Options');
disp('----------------');

% Select Times of Interest
valid = false;
while ~valid
    disp(' ');
    
    selection = input('Restrict Data Range? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        startTime = inputTime('Start');

        if startTime == -1
            continue;
        end

        endTime = inputTime('End');

        if endTime == -1
            continue;
        elseif endTime < startTime
            disp('        WARNING: Invalid Time Format (''endTime'' Precedes ''startTime'')');
            
            continue;
        elseif endTime < expSprayData.time(1) || startTime > expSprayData.time(end)
            disp('        WARNING: No Lagrangian Data in Selected Time Range');
            
            continue;
        end

        i = 1;
        while i <= height(expSprayData.time)

            if expSprayData.time(i) < startTime || expSprayData.time(i) > endTime
                expSprayData.time(i) = [];
                expSprayData.density(i) = [];
            else
                i = i + 1;
            end

        end
        clear i;

        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

% Specify Sampling Frequency
valid = false;
while ~valid
    disp(' ');
    
    disp(['Default Sampling Frequency: ', num2str(sampleFreq), ' Hz']);
    selection = input('    Reduce Recording Frequency? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        sampleInt = 1;

        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        sampleInt = inputFreq(sampleFreq);

        if sampleInt == -1
            continue;
        end

        if sampleInt >= (floor(height(expSprayData.time) / 2))
            disp('            WARNING: Sampling Interval Must Fall Within Data Range');
        else
            valid = true;
        end

    else
        disp('        WARNING: Invalid Entry');
    end

end
clear valid;

if sampleInt ~= 1
    expSprayData.time = expSprayData.time(1:sampleInt:end);
    expSprayData.density = expSprayData.density(1:sampleInt:end);
end

% Define Data ID
timePrecision = 3;

startInst = erase(num2str(expSprayData.time(1), ['%.', num2str(timePrecision), 'f']), '.');
endInst = erase(num2str(expSprayData.time(end), ['%.', num2str(timePrecision), 'f']), '.');

sampleFreq = num2str(sampleFreq / sampleInt);

% Generate Dataset ID
dataID = ['T', startInst, '_T', endInst, '_F', sampleFreq];

% Identify Plane Position
planePos = str2double(planeID(1:(end-1)));

planePos = planePos * normLength;

disp(' ');
disp(' ');


%% Generate Contaminant Maps

disp('Contaminant Mapping');
disp('--------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

%%%%

nTimes = height(expSprayData.time);

disp(' ');

% Map Raw Data Onto Uniform Grid
disp('    Mapping Raw Data Onto Uniform Grid...');

gridShape = [height(unique(expSprayData.positionGrid(:,2))), ...
             height(unique(expSprayData.positionGrid(:,3)))];

xLimsData = planePos;
yLimsData = [-0.43; 0.25];
zLimsData = [0.01; 0.4995];

cellSize.DaVis = double(min(abs(unique(diff(expSprayData.positionGrid(:,2))))));

if strcmp(campaignID, 'Far_Field_Soiling_07_22')
    cellSize.target = 0.8e-3;
end

% Adjust Uniform Cell Size to Fit Region of Interest
nPy = (diff(yLimsData) / cellSize.target) + 1;
nPz = (diff(zLimsData) / cellSize.target) + 1;

sizeY = diff(linspace(yLimsData(1), yLimsData(2), nPy));
sizeZ = diff(linspace(zLimsData(1), zLimsData(2), nPz));

cellSize.y = sizeY(1); clear sizeY;
cellSize.z = sizeZ(1); clear sizeZ;
cellSize.area = cellSize.y * cellSize.z;

% Generate Grid
[y, z] = ndgrid(linspace(yLimsData(1), yLimsData(2), nPy), ...
                linspace(zLimsData(1), zLimsData(2), nPz));

mapData.positionGrid = zeros([height(y(:)),3]);
mapData.positionGrid(:,1) = xLimsData;
mapData.positionGrid(:,(2:3)) = [y(:), z(:)];

nCells = height(mapData.positionGrid);

mapData.time = expSprayData.time;

% Initialise Progress Bar
wB = waitbar(0, 'Mapping Raw Data Onto Uniform Grid', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Perform Mapping Operation
mapData.density.inst = cell(nTimes,1);

densityOrig = expSprayData.density;
yOrig = reshape(expSprayData.positionGrid(:,2), gridShape);
zOrig = reshape(expSprayData.positionGrid(:,3), gridShape);
clear expSprayData
for i = 1:nTimes
    density = reshape(densityOrig{i}, gridShape);
    
    interp = griddedInterpolant(yOrig, zOrig, density, 'linear', 'none');
    
    mapData.density.inst{i} = interp(y, z);
    mapData.density.inst{i} = mapData.density.inst{i}(:);
    
    mapData.density.inst{i}(mapData.density.inst{i} < 1e-6) = 0;
    
    % Update Waitbar
    waitbar((i / nTimes), wB);
end
clear i densityOrig yOrig zOrig y z density;

delete(wB);

disp(' ');

% Calculate Instantaneous Centre of Spray
disp('    Calculating Instantaneous Centre of Spray...');

mapData.CoM.inst = cell(nTimes,1); mapData.CoM.inst(:) = {zeros([1,3], 'single')};

for i = 1:nTimes
    mapData.CoM.inst{i}(1) = mapData.positionGrid(1,1);
    mapData.CoM.inst{i}(2) = sum(mapData.density.inst{i} .* mapData.positionGrid(:,2)) / ...
                                 sum(mapData.density.inst{i});
    mapData.CoM.inst{i}(3) = sum(mapData.density.inst{i} .* mapData.positionGrid(:,3)) / ...
                                 sum(mapData.density.inst{i});
end
clear i;

disp(' ');

% Calculate Time-Averaged Spray Density
disp('    Calculating Time-Averaged Spray Density...');

mapData.density.mean = zeros([nCells,1], 'single');

for i = 1:nTimes
    mapData.density.mean = mapData.density.mean + mapData.density.inst{i};
end
clear i;

mapData.density.mean = mapData.density.mean / nTimes;

disp(' ');

% Calculate Time-Averaged Centre of Spray
disp('    Calculating Time-Averaged Centre of Spray...');

mapData.CoM.mean = zeros([1,3], 'single');

mapData.CoM.mean(1) = mapData.positionGrid(1,1);
mapData.CoM.mean(2) = sum(mapData.density.mean .* mapData.positionGrid(:,2)) / ...
                           sum(mapData.density.mean);
mapData.CoM.mean(3) = sum(mapData.density.mean .* mapData.positionGrid(:,3)) / ...
                          sum(mapData.density.mean);

disp(' ');

% Calculate Instantaneous Spray Density Fluctuations
disp('    Calculating Instantaneous Density Fluctuations...');

mapData.density.prime = mapData.density.inst;

for i = 1:nTimes
    mapData.density.prime{i} = mapData.density.prime{i} - mapData.density.mean;
end
clear i;

disp(' ');

% Calculate RMS of Spray Density
disp('    Calculating RMS of Spray Density...');

mapData.density.RMS = zeros([nCells,1], 'single');

for i = 1:nTimes
    mapData.density.RMS = mapData.density.RMS + mapData.density.prime{i}.^2;
end
clear i;

mapData.density.RMS = sqrt((1 / nTimes) * mapData.density.RMS);

%%%%

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
        
        if ~exist([saveLoc, '/Experimental/MATLAB/planarSprayMap/', campaignID, '/', caseID], 'dir')
            mkdir([saveLoc, '/Experimental/MATLAB/planarSprayMap/', campaignID, '/', caseID]);
        end
        
        disp(['    Saving to: ', saveLoc, '/Experimental/MATLAB/planarSprayMap/', campaignID, '/', caseID, '/', dataID '.mat']);
        save([saveLoc, '/Experimental/MATLAB/planarSprayMap/', campaignID, '/', caseID, '/', dataID '.mat'], ...
             'campaignID', 'caseID', 'planeID', 'dataID', 'mapData', 'cellSize', 'sampleInt', 'timePrecision', 'normDims', '-v7.3', '-noCompression');
        disp('        Success');
        
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
    
    selection = input('Plot Time-Averaged Map? [y/n]: ', 's');

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
    
    selection = input('Plot RMS Map? [y/n]: ', 's');

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
    
    selection = input('Plot Instantaneous Maps? [y/n]: ', 's');

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

        xLimsData = xLimsData / normLength;
        yLimsData = yLimsData / normLength;
        zLimsData = zLimsData / normLength;

        cellSize.DaVis = cellSize.DaVis / normLength;
        cellSize.target = cellSize.target / normLength;
        cellSize.y = cellSize.y / normLength;
        cellSize.z = cellSize.z / normLength;
        cellSize.area = cellSize.area / (normLength^2);

        mapData.positionGrid(:,(2:3)) = mapData.positionGrid(:,(2:3)) / normLength;

        mapData.CoM.mean = mapData.CoM.mean / normLength;

        for i = 1:nTimes
            mapData.CoM.inst{i} = mapData.CoM.inst{i} / normLength;

            % Update Waitbar
            waitbar((i / nTimes), wB);
        end
        clear i;

        delete(wB);
    end
    
    % Normalise Contaminant Maps
    if normDensity
        disp(' ');

        disp('    Normalising Spray Density...');

        wB = waitbar(0, 'Normalising Spray Density', 'name', 'Progress');
        wB.Children.Title.Interpreter = 'none';

        mapData.density.mean = mapData.density.mean / normValue;
        mapData.density.RMS = mapData.density.RMS / normValue;

        for i = 1:nTimes
            mapData.density.inst{i} = mapData.density.inst{i} / normValue;
            mapData.density.prime{i} = mapData.density.prime{i} / normValue;

            % Update Waitbar
            waitbar((i / nTimes), wB);
        end
        clear i;

        delete(wB);
    end
    
end

disp(' ');
disp(' ');


%% Present Contaminant Maps

disp('Map Presentation');
disp('-----------------');

disp(' ');

if plotMean || plotRMS || plotInst
    orientation = 'YZ';
    
    if strcmp(campaignID, 'Far_Field_Soiling_07_22')
        spatialRes = 0.5e-3;
    end
    
    if normDims
        spatialRes = spatialRes / normLength;
    end
    
    positionData = mapData.positionGrid;
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
    
    % Remove Geometry From Empty Tunnel
    if contains(caseID, 'ET')
        geometry = [];
    end
    
end

if plotMean
    disp('    Presenting Time-Averaged Seeding Density...');

    scalarData = mapData.density.mean;
    figName = ['Average_', caseID];
    contourlines = [0.02; 0.02];
    figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
    cLims = [0; 1.05];

    [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                           xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                           planeNo, fig, figName, cMap, geometry, contourlines, ...
                                           refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                           zLimsPlot, normDims, figSave);
    
    disp(' ');
end

if plotRMS
    disp('    Presenting RMS of Seeding Density...');
    
    scalarData = mapData.density.RMS / mean(mapData.density.mean);
    figName = ['RMS_', caseID];
    contourlines = [];
    figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
    cLims = [0; 5];

    [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                           xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                           planeNo, fig, figName, cMap, geometry, contourlines, ...
                                           refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                           zLimsPlot, normDims, figSave);
    
    disp(' ');
end

if plotInst
    disp('    Presenting Instantaneous Seeding Density...');
    
    contourlines = [];
    cLims = [0; 3.6];

    figHold = fig;

    for i = startFrame:endFrame

        if i ~= startFrame
            clf(fig);
            fig = figHold;
        end

        scalarData = mapData.density.inst{i};
        figTime = num2str(mapData.time(i), '%.3f');
        figName = ['Inst_T', erase(figTime, '.'), '_', caseID];
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
    disp('Skipping Map Presentation');
end


%% Local Functions

function time = inputTime(type)

    time = str2double(input(['    Input ', type, ' Time [s]: '], 's'));
    
    if isnan(time) || length(time) > 1 || time <= 0
        disp('        WARNING: Invalid Entry');
        
        time = -1;
    end

end


function sampleInterval = inputFreq(origFreq)
    
    newFreq = str2double(input('        Input Frequency [Hz]: ', 's'));
    
    if isnan(newFreq) || newFreq <= 0 || newFreq > origFreq
        disp('            WARNING: Invalid Entry');
        
        sampleInterval = -1;
    elseif mod(origFreq, newFreq) ~= 0
        disp(['            WARNING: New Frequency Must Be a Factor of ', num2str(origFreq),' Hz']);
        
        sampleInterval = -1;
    else
        sampleInterval = origFreq / newFreq;
    end
    
end


function frameNo = inputFrames(Nt, type)

    frameNo = str2double(input(['    Input Desired ', type, ' Frame [1-', num2str(Nt), ']: '], 's'));
    
    if isnan(frameNo) || frameNo < 1 || frameNo > Nt
        disp('        WARNING: Invalid Entry');
        
        frameNo = -1;
    end

end