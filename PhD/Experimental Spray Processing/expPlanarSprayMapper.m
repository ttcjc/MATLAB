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

if exist('/media/ttcjc/B418689', 'dir')
    dataLocation = '/media/ttcjc/B418689/Data';
else
    dataLocation = '/mnt/Processing/Data';
end

nProc = 4; % Number of Processors Used for Process-Based Parallelisation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures


%% Planar Experimental Spray Mapper v2.0

figSave = false; % Save .fig File(s);

normDims = true; % Normalise Spatial Dimensions

normDensity = true; % Normalise Spray Density
    normValue = 0.0052140; % (SB_1.0L_120s_15Hz_02)

disp('=====================================');
disp('Planar Experimental Spray Mapper v2.0');
disp('=====================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Update to Support Changes to 'plotPlanarScalarField'
% v2.0 - Update To Support Changes to DaVis Data Format


%% Initialise Case

% Select Relevant Geometry and Define Bounding Box
[geometry, xDims, yDims, zDims, spacePrecision, normDims, normLength] = selectGeometry(normDims);

disp(' ');
disp(' ');


%% Initialise Experimental Spray Data

[campaignID, caseID, planeID, ...
 expSprayData, sampleFreq] = initialiseExpSprayData(saveLocation, dataLocation, nProc);

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

if normDims
    dataID = ['T', startInst, '_T', endInst, '_F', sampleFreq, '_Norm'];
else
    dataID = ['T', startInst, '_T', endInst, '_F', sampleFreq];
end

% Identify Plane Position
planePos = str2double(planeID(1:(end-1)));

if ~normDims
    
    if strcmp(campaignID, 'Far_Field_Soiling_07_22')
        planePos = planePos * 1.044;
    end
    
end

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

% Map Raw Data Onto Adjusted Grid
disp('    Mapping Raw Data Onto Adjusted Grid...');

gridShape = [height(unique(expSprayData.positionGrid(:,2))), ...
             height(unique(expSprayData.positionGrid(:,3)))];

xLimsData = planePos;
yLimsData = [-0.43; 0.25];
zLimsData = [0.01; 0.4995];

cellSize.DaVis = double(min(abs(unique(diff(expSprayData.positionGrid(:,2))))));

if strcmp(campaignID, 'Far_Field_Soiling_07_22')
    cellSize.target = 0.8e-3;
end

cellSize.y = (yLimsData (2) - yLimsData (1)) / round(((yLimsData (2) - yLimsData (1)) / cellSize.target));
cellSize.z = (zLimsData (2) - zLimsData (1)) / round(((zLimsData (2) - zLimsData (1)) / cellSize.target));

[y, z] = ndgrid(yLimsData(1):cellSize.y:yLimsData(2), zLimsData(1):cellSize.z:zLimsData(2));

mapData.positionGrid = zeros([height(y(:)),3]);
mapData.positionGrid(:,1) = xLimsData;
mapData.positionGrid(:,(2:3)) = [y(:), z(:)];

mapData.time = expSprayData.time;

% Initialise Progress Bar
wB = waitbar(0, 'Mapping Raw Data Onto Adjusted Grid', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Perform Mapping Operation
mapData.density.inst = cell(nTimes,1);

yOrig = reshape(expSprayData.positionGrid(:,2), gridShape);
zOrig = reshape(expSprayData.positionGrid(:,3), gridShape);
for i = 1:nTimes
    density = reshape(expSprayData.density{i}, gridShape);
    
    interp = griddedInterpolant(yOrig, zOrig, density, 'linear', 'none');
    
    mapData.density.inst{i} = interp(y, z);
    mapData.density.inst{i} = mapData.density.inst{i}(:);
    
    mapData.density.inst{i}(mapData.density.inst{i} < 1e-6) = 0;
    
    % Update Waitbar
    waitbar((i / nTimes), wB);
end
clear i yOrig zOrig y z density expSprayData;

delete(wB);

% Normalise Position Data
if normDims
    mapData.positionGrid(:,(2:3)) = round((mapData.positionGrid(:,(2:3)) / normLength), spacePrecision);

    yLimsData = round((yLimsData / normLength), spacePrecision);
    zLimsData = round((zLimsData / normLength), spacePrecision);
end

disp(' ');

% Normalise Contaminant Maps
if normDensity
    disp('    Normalising Instantaneous Spray Density...');
    
    for i = 1:nTimes %#ok<*UNRCH>
        mapData.density.inst{i} = mapData.density.inst{i} / normValue;
    end
    clear i;
    
end

disp(' ');

% Calculate Time-Averaged Spray Density
disp('    Calculating Time-Averaged Spray Density...');

mapData.density.mean = zeros([height(mapData.positionGrid),1], 'single');

for i = 1:nTimes
    mapData.density.mean = mapData.density.mean + mapData.density.inst{i};
end
clear i;

mapData.density.mean = mapData.density.mean / nTimes;

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

mapData.density.RMS = zeros([height(mapData.positionGrid),1]);

for i = 1:nTimes
    mapData.density.RMS = mapData.density.RMS + mapData.density.prime{i}.^2;
end
clear i;

mapData.density.RMS = sqrt((1 / nTimes) * mapData.density.RMS);

disp(' ');

% Calculate Instantaneous Centre of Spray
disp('    Calculating Instantaneous Centre of Spray...');

mapData.CoM.inst = cell(nTimes,1); mapData.CoM.inst(:) = {zeros([1,3], 'single')};

for i = 1:nTimes
    mapData.CoM.inst{i}(1) = mapData.positionGrid(1,1);
    mapData.CoM.inst{i}(2) = sum(mapData.density.inst{i} .* mapData.positionGrid(:,2)) / sum(mapData.density.inst{i});
    mapData.CoM.inst{i}(3) = sum(mapData.density.inst{i} .* mapData.positionGrid(:,3)) / sum(mapData.density.inst{i});
end
clear i;

disp(' ');

% Calculate Time-Averaged Centre of Spray
disp('    Calculating Time-Averaged Centre of Spray...');

mapData.CoM.mean = zeros([1,3], 'single');

mapData.CoM.mean(1) = mapData.positionGrid(1,1);
mapData.CoM.mean(2) = sum(mapData.density.mean .* mapData.positionGrid(:,2)) / sum(mapData.density.mean);
mapData.CoM.mean(3) = sum(mapData.density.mean .* mapData.positionGrid(:,3)) / sum(mapData.density.mean);

%%%%

executionTime = toc;

disp(' ');

disp(['    Run Time: ', num2str(executionTime), 's']);

disp(' ');

disp('  SUCCESS  ');
disp('***********');

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

disp(' ');
disp(' ');


%% Present Contaminant Maps

disp('Map Presentation');
disp('-----------------');

disp(' ');

if plotMean || plotRMS || plotInst
    orientation = 'YZ';
    
    if strcmp(campaignID, 'Far_Field_Soiling_07_22')
        xLimsPlot = [0.31875; 2.73575];
        yLimsPlot = [-0.522; 0.522];
        zLimsPlot = [0; 0.522];
    end

    if normDims
        xLimsPlot = round((xLimsPlot / normLength), spacePrecision);
        yLimsPlot = round((yLimsPlot / normLength), spacePrecision);
        zLimsPlot = round((zLimsPlot / normLength), spacePrecision);
    end
    
    if strcmp(campaignID, 'Far_Field_Soiling_07_22')
        spatialRes = 0.5e-3;
    end
    
    positionData = mapData.positionGrid;
    mapPerim = [];
    nPlanes = 1;
    planeNo = 1;
    cMap = flipud(viridis(32));
    refPoint = [];
end

% Remove Geometry From Empty Tunnel
if contains(caseID, 'ET')
    geometry = [];
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
        figName = ['Instantaneous_T', erase(figTime, '.'), '_', caseID];
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
        
        if ~exist([saveLocation, '/Experimental/MATLAB/planarSprayMap/', caseID], 'dir')
            mkdir([saveLocation, '/Experimental/MATLAB/planarSprayMap/', caseID]);
        end
        
        disp(['    Saving to: ', saveLocation, '/Experimental/MATLAB/planarSprayMap/', caseID, '/', dataID '.mat']);
        save([saveLocation, '/Experimental/MATLAB/planarSprayMap/', caseID, '/', dataID '.mat'], ...
             'caseID', 'planeID', 'dataID', 'mapData', 'cellSize', 'sampleInt', 'timePrecision', 'normDims', '-v7.3', '-noCompression');
        disp('        Success');
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;


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
