%% Lagrangian Line of Sight Calculator v1.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

saveLocation = '/mnt/Processing/Data';
% saveLocation = '~/Data';

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

samplesPerCell = 32; % Number of Samples per Cell Used During LOS Interpolation

disp('=============================');
disp('Line of Sight Calculator v1.0');
disp('=============================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit


%% Select Region of Interest

disp('Region of Interest');
disp('-------------------');

disp(' ');

disp('Possible Regions of Interest:');
disp('    A: Near-Field');
disp('    B: Far-Field');

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
    [fileName, filePath] = uigetfile([saveLocation, '/Numerical/MATLAB/volumeField/*.mat'], ...
                                      'Select Volumetric Data');
    
    switch format
        
        case 'A'
            
            if contains(filePath, '/nearField')
                disp(['Loading ''', fileName, '''...']);
                dataID = load([filePath, fileName], 'dataID').dataID;
                volumeData = load([filePath, fileName], 'volumeData').volumeData;
                cellSize = load([filePath, fileName], 'cellSize').cellSize;
                sampleInterval = load([filePath, fileName], 'sampleInterval').sampleInterval;
                dLims = load([filePath, fileName], 'dLims').dLims;
                normalise = load([filePath, fileName], 'normalise').normalise;
                disp('    Success');
                
                valid = true;
            else
                disp('WARNING: Invalid File Selection');
                clear fileName filePath;
            end
            
        case 'B'
            
            if contains(filePath, '/farField')
                disp(['Loading ''', fileName, '''...']);
                dataID = load([filePath, fileName], 'dataID').dataID;
                volumeData = load([filePath, fileName], 'volumeData').volumeData;
                cellSize = load([filePath, fileName], 'cellSize').cellSize;
                sampleInterval = load([filePath, fileName], 'sampleInterval').sampleInterval;
                dLims = load([filePath, fileName], 'dLims').dLims;
                normalise = load([filePath, fileName], 'normalise').normalise;
                disp('    Success');
                
                valid = true;
            else
                disp('WARNING: Invalid File Selection');
                clear fileName filePath;
            end
            
    end
end
clear valid;

namePos = strfind(filePath, '/');
caseName = filePath((namePos(end - 2) + 1):(namePos(end - 1) - 1));

timePrecision = strfind(fileName, '_T') - 3;

disp(' ');
disp(' ');


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normalise] = selectGeometry(normalise);

disp(' ');
disp(' ');


%% Select Plane of Interest

% Temporarily Restore Original Dimensions
if normalise

    if contains(caseName, ["Run_Test", "Windsor"])
        volumeData.positionGrid = round((volumeData.positionGrid * 1.044), spacePrecision);
    end

end

% Select Plane of Interest
LOSdata = identifyVolumeSlices(volumeData.positionGrid, spacePrecision, false);

% Extract Planar Position Data
index = find(volumeData.positionGrid(:,1) == LOSdata.position);
LOSdata.positionGrid = volumeData.positionGrid(index,:);

disp(' ');
disp(' ');


%% Select Origin Point

disp('Origin Point Definition');
disp('------------------------');

% Select Origin Point
LOSdata.originPoint = zeros(1,3);

valid = false;
while ~valid
    disp(' ')
    disp('Specify Origin Point:')
    
    LOSdata.originPoint(1) = inputPos('X');
    LOSdata.originPoint(2) = inputPos('Y');
    LOSdata.originPoint(3) = inputPos('Z');

    if (LOSdata.originPoint(1) < min(volumeData.positionGrid(:,1)) || LOSdata.originPoint(1) > max(volumeData.positionGrid(:,1))) || ...
       (LOSdata.originPoint(2) < min(volumeData.positionGrid(:,2)) || LOSdata.originPoint(2) > max(volumeData.positionGrid(:,1))) || ...
       (LOSdata.originPoint(3) < min(volumeData.positionGrid(:,3)) || LOSdata.originPoint(3) > max(volumeData.positionGrid(:,1)))
        disp('        WARNING: Origin Point Lies Outside Volume');
        continue;
    end
    
    switch LOSdata.orientation
        
        case 'YZ'
            
            if LOSdata.originPoint(1) < LOSdata.position
                disp('        WARNING: Origin Point Must Lie in the Positive X-Direction of the Target Plane');
                continue;
            end
            
        case 'XZ'
            
            if LOSdata.originPoint(2) > LOSdata.position
                disp('        WARNING: Origin Point Must Lie in the Negative Y-Direction of the Target Plane');
                continue;
            end
            
        case 'XY'
            
            if LOSdata.originPoint(3) > LOSdata.position
                disp('        WARNING: Origin Point Must Lie in the Positive Z-Direction of the Target Plane');
                continue;
            end
            
    end
    
    valid = true;    
end
clear valid;

% Restore Normalised Dimensions
if normalise

    if contains(caseName, ["Run_Test", "Windsor"])
        volumeData.positionGrid = round((volumeData.positionGrid / 1.044), spacePrecision);
        LOSdata.originPoint = round((LOSdata.originPoint / 1.044), spacePrecision);
        LOSdata.position = round((LOSdata.position / 1.044), spacePrecision);
    end

end

disp(' ');
disp(' ');


%% Calculate Line of Sight

disp('Line of Slight Calculation');
disp('---------------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

disp(' ');

disp('    Initialising...');

evalc('parpool(nProc);');

disp(' ');

% Check if Plane of Interest Intersects Geometry
removeIntersect = false;

switch LOSdata.orientation
    
    case 'YZ'
        
        if LOSdata.position <= xDims(2)
            removeIntersect = true;
        end
        
    case 'XZ'
        
        if LOSdata.position >= yDims(1)
            removeIntersect = true;
        end
        
    case 'XY'
        
        if LOSdata.position <= zDims(2)
            removeIntersect = true;
        end
        
end

if removeIntersect
    disp('    Removing Erroneous Data From Volume Nodes Intersecting Geometry...');

    % Remove Erroneous Data From Cells Intersecting Geometry
    parts = fieldnames(geometry);
    fields = fieldnames(volumeData.mean);
    for i = 1:height(parts)
        geoPoints = unique(geometry.(parts{i}).vertices, 'rows');
        DT = delaunay(geoPoints);

        index = ~isnan(tsearchn(geoPoints, DT, volumeData.positionGrid));

        for j = 1:height(fields)
            volumeData.mean.(fields{j})(index,:) = nan;

            for k = 1:height(volumeData.inst.time)
                volumeData.inst.(fields{j}){k}(index,:) = nan;
            end

        end

    end
    clear parts fields;

    disp(' ');
end

disp('    Reshaping Position Data for Improved Interpolation Performance...');

% Reshape Position Data for Improved Interpolation Performance
gridShape = [height(unique(volumeData.positionGrid(:,1))), ...
             height(unique(volumeData.positionGrid(:,2))), ...
             height(unique(volumeData.positionGrid(:,3)))];

x = reshape(volumeData.positionGrid(:,1), gridShape);
y = reshape(volumeData.positionGrid(:,2), gridShape);
z = reshape(volumeData.positionGrid(:,3), gridShape);

disp(' ');

disp('    Calculating Instantaneous Line of Sight...');

LOSdata.inst.time = volumeData.inst.time;

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Line of Sight', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, height(LOSdata.inst.time));

% Calculate Instantaneous Line of Sight
massInLOS = cell(height(LOSdata.inst.time),1);

mass = volumeData.inst.mass;
cellVolume = cellSize.volume;
positionGrid = LOSdata.positionGrid;
orientation = LOSdata.orientation;
cellSizeX = cellSize.x;
cellSizeY = cellSize.y;
cellSizeZ = cellSize.z;
originPointX = LOSdata.originPoint(1);
originPointY = LOSdata.originPoint(2);
originPointZ = LOSdata.originPoint(3);
position = LOSdata.position;
parfor i = 1:height(LOSdata.inst.time)
    sprayDensity = reshape((mass{i} / cellVolume), gridShape);
    
    interp = griddedInterpolant(x, y, z, sprayDensity, 'linear', 'none');
    
    massInLOS{i} = zeros(height(positionGrid),1);
    
    switch orientation
        
        case 'YZ'
            dX = cellSizeX / samplesPerCell;
            sampleX = (originPointX:-dX:position)';
            
            for j = 1:height(positionGrid)
                t = (sampleX - originPointX) / (positionGrid(j,1) - originPointX);
                sampleXYZ = ([originPointX, originPointY, originPointZ] + t*(positionGrid(j,:) - [originPointX, originPointY, originPointZ]));
                interpMass = interp(sampleXYZ(:,1), sampleXYZ(:,2), sampleXYZ(:,3)) * dX;
                
                if sum(isnan(interpMass)) > samplesPerCell
                    massInLOS{i}(j) = nan;
                else
                    massInLOS{i}(j) = sum(interpMass, 'omitNaN');
                end
                
            end
            
        case 'XZ'
            dY = cellSizeY / samplesPerCell;
            sampleY = (originPointY:dY:position)';
            
            for j = 1:height(positionGrid)
                t = (sampleY -originPointY) / (positionGrid(j,2) - originPointY);
                sampleXYZ = ([originPointX, originPointY, originPointZ] + t*(positionGrid(j,:) - [originPointX, originPointY, originPointZ]));
                interpMass = interp(sampleXYZ(:,1), sampleXYZ(:,2), sampleXYZ(:,3)) * dY;
                
                if sum(isnan(interpMass)) > samplesPerCell
                    massInLOS{i}(j) = nan;
                else
                    massInLOS{i}(j) = sum(interpMass, 'omitNaN');
                end
                
            end
            
        case 'XY'
            dZ = cellSizeZ / samplesPerCell;
            sampleZ = (originPointZ:-dZ:position)';
            
            for j = 1:height(positionGrid)
                t = (sampleZ - originPointZ) / (positionGrid(j,3) - originPointZ);
                sampleXYZ = ([originPointX, originPointY, originPointZ] + t*(positionGrid(j,:) - [originPointX, originPointY, originPointZ]));
                interpMass = interp(sampleXYZ(:,1), sampleXYZ(:,2), sampleXYZ(:,3)) * dZ;
                
                if sum(isnan(interpMass)) > samplesPerCell
                    massInLOS{i}(j) = nan;
                else
                    massInLOS{i}(j) = sum(interpMass, 'omitNaN');
                end
                
            end
            
    end
    
    send(dQ, []);
end
clear mass cellVolume positionGrid orientation cellSizeX cellSizeY cellSizeZ originPointX originPointY originPointZ position;

delete(wB);

LOSdata.inst.mass = massInLOS;

clear massInLOS;

disp(' ');

disp('    Calculating Time-Averaged Line of Sight...');

% Calculate Time-Averaged Line of Sight
sprayDensity = reshape((volumeData.mean.mass / cellSize.volume), gridShape);

interp = griddedInterpolant(x, y, z, sprayDensity, 'linear', 'none');

LOSdata.mean.mass = zeros(height(LOSdata.positionGrid),1);

switch LOSdata.orientation

    case 'YZ'
        dX = cellSize.x / samplesPerCell;
        sampleX = (LOSdata.originPoint(1):-dX:LOSdata.position)';

        for i = 1:height(LOSdata.positionGrid)
            t = (sampleX - LOSdata.originPoint(1)) / (LOSdata.positionGrid(i,1) - LOSdata.originPoint(1));
            sampleXYZ = (LOSdata.originPoint + t*(LOSdata.positionGrid(i,:) - LOSdata.originPoint));
            interpMass = interp(sampleXYZ(:,1), sampleXYZ(:,2), sampleXYZ(:,3)) * dX;
            
            if sum(isnan(interpMass)) > samplesPerCell
                LOSdata.mean.mass(i) = nan;
            else
                LOSdata.mean.mass(i) = sum(interpMass, 'omitNaN');
            end
            
        end

    case 'XZ'
        dY = cellSize.y / samplesPerCell;
        sampleY = (LOSdata.originPoint(2):dY:LOSdata.position)';

        for i = 1:height(LOSdata.positionGrid)
            t = (sampleY - LOSdata.originPoint(2)) / (LOSdata.positionGrid(i,2) - LOSdata.originPoint(2));
            sampleXYZ = (LOSdata.originPoint + t*(LOSdata.positionGrid(i,:) - LOSdata.originPoint));
            interpMass = interp(sampleXYZ(:,1), sampleXYZ(:,2), sampleXYZ(:,3)) * dY;
            
            if sum(isnan(interpMass)) > samplesPerCell
                LOSdata.mean.mass(i) = nan;
            else
                LOSdata.mean.mass(i) = sum(interpMass, 'omitNaN');
            end
            
        end

    case 'XY'
        dZ = cellSize.z / samplesPerCell;
        sampleZ = (LOSdata.originPoint(3):-dZ:LOSdata.position)';

        for i = 1:height(LOSdata.positionGrid)
            t = (sampleZ - LOSdata.originPoint(3)) / (LOSdata.positionGrid(i,3) - LOSdata.originPoint(3));
            sampleXYZ = (LOSdata.originPoint + t*(LOSdata.positionGrid(i,:) - LOSdata.originPoint));
            interpMass = interp(sampleXYZ(:,1), sampleXYZ(:,2), sampleXYZ(:,3)) * dZ;
            
            if sum(isnan(interpMass)) > samplesPerCell
                LOSdata.mean.mass(i) = nan;
            else
                LOSdata.mean.mass(i) = sum(interpMass, 'omitNaN');
            end
            
            
        end

end

evalc('delete(gcp(''nocreate''));');

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
    selection = input('Plot Time-Averaged Map(s)? [y/n]: ', 's');

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
    selection = input('Plot Instantaneous Maps? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotInst = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotInst = true;
        nFrames = inputFrames(height(LOSdata.inst.time));
        
        if nFrames == -1
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

if plotInst || plotMean
    
    switch LOSdata.orientation

        case 'YZ'
            xLimsData = LOSdata.position;
            yLimsData = [min(LOSdata.positionGrid(:,2)); max(LOSdata.positionGrid(:,2))];
            zLimsData = [min(LOSdata.positionGrid(:,3)); max(LOSdata.positionGrid(:,3))];

        case 'XZ'
            xLimsData = [min(LOSdata.positionGrid(:,1)); max(LOSdata.positionGrid(:,1))];
            yLimsData = LOSdata.position;
            zLimsData = [min(LOSdata.positionGrid(:,3)); max(LOSdata.positionGrid(:,3))];

        case 'XY'
            xLimsData = [min(LOSdata.positionGrid(:,1)); max(LOSdata.positionGrid(:,1))];
            yLimsData = [min(LOSdata.positionGrid(:,2)); max(LOSdata.positionGrid(:,2))];

    end
    
    switch format

        case 'A'

            if contains(caseName, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 1.26625]; % 0.75 L
                yLimsPlot = [-0.3445; 0.3445];
                zLimsPlot = [0; 0.489];
            end

        case 'B'

            if contains(caseName, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 2.57125]; % 2 L
                yLimsPlot = [-0.5945; 0.5945];
                zLimsPlot = [0; 0.739];
%                 xLimsPlot = [0.31875; 2.57125]; % 2 L
%                 yLimsPlot = [-0.4945; 0.4945];
%                 zLimsPlot = [0; 0.639];
            end

    end

    if normalise
        xLimsPlot = round((xLimsPlot / 1.044), spacePrecision);
        yLimsPlot = round((yLimsPlot / 1.044), spacePrecision);
        zLimsPlot = round((zLimsPlot / 1.044), spacePrecision);
    end
    
    orientation = LOSdata.orientation;
    positionData = LOSdata.positionGrid;
    mapPerim = [];
    cMap = flipud(viridis(32));
    contourlines = [];
    CoM = LOSdata.originPoint;
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    nPlanes = 1;
    planeNo = 1;
    
end

if plotMean
    disp('    Presenting Time-Averaged Line of Sight Map...');
    
    scalarData = LOSdata.mean.mass;
    figName = 'Time_Averaged_LOS_Map';
    figSubtitle = ' ';
%     cLims = [0; max(LOSdata.mean.mass)];
    cLims = [0; 5e-3];
    
    [fig, planeNo] = planarScalarPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                       mapPerim, fig, figName, cMap, geometry, contourlines, ...
                                       xDims, yDims, zDims, CoM, figTitle, figSubtitle, cLims, ...
                                       xLimsPlot, yLimsPlot, zLimsPlot, normalise, nPlanes, planeNo);
    
    disp(' ');
end

if plotInst
    disp('    Presenting Instantaneous Line of Sight Map(s)...');
    
%     cLims = [0; max(cellfun(@max, LOSdata.inst.mass))];
    cLims = [0; 5e-3];
    
    figHold = fig;
    
    for i = 1:nFrames
        
        if i ~= 1
            clf(fig);
            fig = figHold;
        end
        
        scalarData = LOSdata.inst.mass{i};
        figTime = num2str(LOSdata.inst.time(i), ['%.', num2str(timePrecision), 'f']);
        figName = ['Instantaneous_LOS_Map_T_', erase(figTime, '.')];
        figSubtitle = [figTime, ' \it{s}'];

        [fig, planeNo] = planarScalarPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                           mapPerim, fig, figName, cMap, geometry, contourlines, ...
                                           xDims, yDims, zDims, CoM, figTitle, figSubtitle, cLims, ...
                                           xLimsPlot, yLimsPlot, zLimsPlot, normalise, nPlanes, planeNo);
    end
    
    disp(' ');
end

if ~plotMean && ~plotInst
    disp('    Skipping Line of Sight Presentation');
    disp(' ');
end

disp(' ');


%% Save Line of Sight Data




%% Local Functions

function pos = inputPos(orientation)

    pos = str2double(input(['    ', orientation, '-Position [m]: '], 's'));
    
    if isnan(pos) || length(pos) > 1
        disp('        WARNING: Invalid Entry');
        pos = -1;
    end
    
end


function nFrames = inputFrames(Nt)

    nFrames = str2double(input(['    Input Desired Frame Count [1-', num2str(Nt), ']: '], 's'));
    
    if isnan(nFrames) || nFrames <= 0 || nFrames > Nt
        disp('        WARNING: Invalid Entry');
        nFrames = -1;
    end

end