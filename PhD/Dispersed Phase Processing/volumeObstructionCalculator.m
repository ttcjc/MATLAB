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


%% Lagrangian Line of Sight Obstruction Calculator v2.0

dL = 1.25e-4; % Small Distance Used During LOS Numerical Integration [m]

figSave = false; % Save .fig File(s)

disp('=========================================');
disp('Line of Sight Obstruction Calculator v2.0');
disp('=========================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v2.0 - 


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
                
                caseID = load([filePath, fileName], 'caseID').caseID;
                dataID = load([filePath, fileName], 'dataID').dataID;
                volumeData = load([filePath, fileName], 'volumeData').volumeData;
                cellSize = load([filePath, fileName], 'cellSize').cellSize;
                sampleInterval = load([filePath, fileName], 'sampleInterval').sampleInterval;
                timePrecision = load([filePath, fileName], 'timePrecision').timePrecision;
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
                
                caseID = load([filePath, fileName], 'caseID').caseID;
                dataID = load([filePath, fileName], 'dataID').dataID;
                volumeData = load([filePath, fileName], 'volumeData').volumeData;
                cellSize = load([filePath, fileName], 'cellSize').cellSize;
                sampleInterval = load([filePath, fileName], 'sampleInterval').sampleInterval;
                timePrecision = load([filePath, fileName], 'timePrecision').timePrecision;
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

disp(' ');
disp(' ');


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normalise] = selectGeometry(normalise);

disp(' ');
disp(' ');


%% Select Plane of Interest

% Temporarily Restore Original Dimensions
if normalise

    if contains(caseID, ["Run_Test", "Windsor"])
        volumeData.positionGrid = round((volumeData.positionGrid * 1.044), spacePrecision);
    end

end

% Select Plane of Interest
obstructData.targetPlane = identifyVolumeSlices(volumeData.positionGrid, spacePrecision, false);

% Extract Planar Position Data
index = find(volumeData.positionGrid(:,1) == obstructData.targetPlane.position);
obstructData.positionGrid = volumeData.positionGrid(index,:);

disp(' ');
disp(' ');


%% Select Origin Point

disp('Ray Origin Definition');
disp('----------------------');

% Select Ray Origin Point
obstructData.originPoint = zeros([1,3]);

valid = false;
while ~valid
    disp(' ')
    disp('Specify Ray Origin Point:')
    
    obstructData.originPoint(1) = inputPos('X');
    obstructData.originPoint(2) = inputPos('Y');
    obstructData.originPoint(3) = inputPos('Z');

    if (obstructData.originPoint(1) < min(volumeData.positionGrid(:,1)) || obstructData.originPoint(1) > max(volumeData.positionGrid(:,1))) || ...
       (obstructData.originPoint(2) < min(volumeData.positionGrid(:,2)) || obstructData.originPoint(2) > max(volumeData.positionGrid(:,1))) || ...
       (obstructData.originPoint(3) < min(volumeData.positionGrid(:,3)) || obstructData.originPoint(3) > max(volumeData.positionGrid(:,1)))
        disp('        WARNING: Origin Point Lies Outside Volume');
        continue;
    end
    
    switch obstructData.targetPlane.orientation
        
        case 'YZ'
            
            if obstructData.originPoint(1) < obstructData.targetPlane.position
                disp('        WARNING: Origin Point Must Lie in the Positive X-Direction of the Target Plane');
                continue;
            end
            
        case 'XZ'
            
            if obstructData.originPoint(2) > obstructData.targetPlane.position
                disp('        WARNING: Origin Point Must Lie in the Negative Y-Direction of the Target Plane');
                continue;
            end
            
        case 'XY'
            
            if obstructData.originPoint(3) > obstructData.targetPlane.position
                disp('        WARNING: Origin Point Must Lie in the Positive Z-Direction of the Target Plane');
                continue;
            end
            
    end
    
    valid = true;    
end
clear valid;

% Restore Normalised Dimensions
if normalise

    if contains(caseID, ["Run_Test", "Windsor"])
        volumeData.positionGrid = round((volumeData.positionGrid / 1.044), spacePrecision);
        obstructData.positionGrid = round((obstructData.positionGrid / 1.044), spacePrecision);
        obstructData.targetPlane.position = round((obstructData.targetPlane.position / 1.044), spacePrecision);
        obstructData.originPoint = round((obstructData.originPoint / 1.044), spacePrecision);
        
        dL = round((dL / 1.044), spacePrecision);
    end

end

disp(' ');
disp(' ');


%% Calculate Line of Sight Obstruction

disp('Line of Slight Obstruction Calculation');
disp('---------------------------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;
evalc('parpool(''threads'');');

%%%%

disp(' ');

disp('    Initialising...');

nTimes = height(volumeData.time);

% Check if Plane of Interest Intersects Geometry
removeIntersect = false;

switch obstructData.targetPlane.orientation
    
    case 'YZ'
        
        if obstructData.targetPlane.position <= xDims(2)
            removeIntersect = true;
        end
        
    case 'XZ'
        
        if obstructData.targetPlane.position >= yDims(1)
            removeIntersect = true;
        end
        
    case 'XY'
        
        if obstructData.targetPlane.position <= zDims(2)
            removeIntersect = true;
        end
        
end

% Remove Erroneous Data From Cells Intersecting Geometry
if removeIntersect
    disp('        Removing Erroneous Data From Grid Cells Intersecting Geometry...');

    % Perform Removal
    parts = fieldnames(geometry);
    fields = fieldnames(volumeData.mean);
    for i = 1:height(parts)
        geoPoints = unique(geometry.(parts{i}).vertices, 'rows');
        DT = delaunay(geoPoints);

        index = ~isnan(tsearchn(geoPoints, DT, volumeData.positionGrid));

        for j = 1:height(fields)
            volumeData.mean.(fields{j})(index,:) = nan;

            for k = 1:nTimes
                volumeData.inst.(fields{j}){k}(index,:) = nan;
            end

        end
        clear j;

    end
    clear i;
    
end

% Identify Rays of Interest for Lidar Calculations
obstructData.raysOfInterest = {'lowerLeft'; 'upperLeft'; 'lowerRight'; 'upperRight'; 'Central'};

rayLL = find((obstructData.positionGrid(:,2) == min(obstructData.positionGrid(:,2)) & ...
                 (obstructData.positionGrid(:,3) == min(obstructData.positionGrid(:,3)))));
rayUL = find((obstructData.positionGrid(:,2) == min(obstructData.positionGrid(:,2)) & ...
                 (obstructData.positionGrid(:,3) == max(obstructData.positionGrid(:,3)))));
rayLR = find((obstructData.positionGrid(:,2) == max(obstructData.positionGrid(:,2)) & ...
                 (obstructData.positionGrid(:,3) == min(obstructData.positionGrid(:,3)))));
rayUR = find((obstructData.positionGrid(:,2) == max(obstructData.positionGrid(:,2)) & ...
                 (obstructData.positionGrid(:,3) == max(obstructData.positionGrid(:,3)))));

if contains(caseID, ["Run_Test", "Windsor"])
    rayCent = dsearchn(obstructData.positionGrid(:,[2,3]), [0, 0.1945]);
    rayInj = dsearchn(obstructData.positionGrid(:,[2,3]), [-0.167 0.006]);
end

[rayID, index] = sort([rayLL; rayUL; rayLR; rayUR; rayCent]); clear rayLL rayUL rayLR rayUR rayCent;
obstructData.raysOfInterest = obstructData.raysOfInterest(index);

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
pointsRay = cell(nTimes,1); pointsRay(:) = {cell(height(rayID),1)};
massRay = pointsRay;
nParticlesRay = pointsRay;
dRay = pointsRay;
areaDensity = cell(nTimes,1); areaDensity(:) = {zeros([height(obstructData.positionGrid),1], 'single')};

mass = volumeData.inst.mass;
nParticles = volumeData.inst.nParticles;
d = volumeData.inst.d20;
positionGrid = obstructData.positionGrid;
originPoint = obstructData.originPoint;
cellVolume = cellSize.volume;
parfor i = 1:nTimes
    massField = reshape(mass{i}, gridShape);
    nParticlesField = reshape(nParticles{i}, gridShape);
    dField = reshape(d{i}, gridShape);
    
    massInterp = griddedInterpolant(x, y, z, massField, 'linear', 'none');
    nParticlesInterp = griddedInterpolant(x, y, z, nParticlesField, 'linear', 'none');
    dInterp = griddedInterpolant(x, y, z, dField, 'nearest', 'none');
    
    k = 1;
    for j = 1:height(positionGrid)
        dirVec = positionGrid(j,:) - originPoint;
        distFull = sqrt(dirVec(1)^2 + dirVec(2)^2 + dirVec(3)^2);

        dist = (dL:dL:distFull)';
        samplePoints = originPoint + (dist * (dirVec / distFull));
        
        areaDensity{i}(j) = sum((massInterp(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3))) * ...
                                (dL / cellVolume));
        
        if j == rayID(k) %#ok<PFBNS>
            pointsRay{i}{k} = samplePoints;
            massRay{i}{k} = massInterp(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3));
            nParticlesRay{i}{k} = nParticlesInterp(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3));
            dRay{i}{k} = dInterp(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3));
            
            if j ~= rayID(end)
                k = k + 1;
            end
            
        end
        
    end
    
    % Remove Unnecessary Data
    mass{i} = [];
    nParticle{i} = [];
    d{i} = [];
    
    % Update Waitbar
    send(dQ, []);
end
clear i j k mass nParticles d positionGrid originPoint cellVolume;

delete(wB);

obstructData.inst.rayData.samplePoints = pointsRay; clear pointsRay;
obstructData.inst.rayData.mass = massRay; clear massRay;
obstructData.inst.rayData.nParticles = nParticlesRay; clear nParticlesRay;
obstructData.inst.rayData.d = dRay; clear dRay;
obstructData.inst.areaDensity = areaDensity; clear areaDensity;

disp(' ');

disp('    Calculating Time-Averaged Obstruction...');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Time-Averaged Obstruction', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Perform Calculation
obstructData.mean.rayData.samplePoints = cell(height(rayID),1);
obstructData.mean.rayData.mass = obstructData.mean.rayData.samplePoints;
obstructData.mean.rayData.nParticles = obstructData.mean.rayData.samplePoints;
obstructData.mean.rayData.d = obstructData.mean.rayData.samplePoints;
obstructData.mean.areaDensity = zeros([height(obstructData.positionGrid),1], 'single');

massField = reshape(volumeData.mean.mass, gridShape);
nParticlesField = reshape(volumeData.mean.nParticles, gridShape);
dField = reshape(volumeData.mean.d20, gridShape);
    
massInterp = griddedInterpolant(x, y, z, massField, 'linear', 'none');
nParticlesInterp = griddedInterpolant(x, y, z, nParticlesField, 'linear', 'none');
dInterp = griddedInterpolant(x, y, z, dField, 'linear', 'none');

j = 1;
for i = 1:height(obstructData.positionGrid)
    dirVec = obstructData.positionGrid(i,:) - obstructData.originPoint;
    distFull = sqrt(dirVec(1)^2 + dirVec(2)^2 + dirVec(3)^2);

    dist = (dL:dL:distFull)';
    samplePoints = obstructData.originPoint + (dist * (dirVec / distFull));
    
    obstructData.mean.areaDensity(i) = sum((massInterp(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3))) * ...
                                      (dL / cellSize.volume));
    
    if i == rayID(j)
        obstructData.mean.rayData.samplePoints{j} = samplePoints;
        obstructData.mean.rayData.mass{j} = massInterp(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3));
        obstructData.mean.rayData.nParticles{j} = nParticlesInterp(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3));
        obstructData.mean.rayData.d{j} = dInterp(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3));
        
        if i ~= rayID(end)
            j = j + 1;
        end
        
    end
    
    % Update Waitbar
    waitbar((i / height(obstructData.positionGrid)), wB);
end
clear i j;

delete(wB);

clear volumeData;

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


%% Select Presentation Options

disp('Presentation Options');
disp('---------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Plot Ray Perspective? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotPerspective = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotPerspective = true;
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

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
    selection = input('Plot Instantaneous Obstruction? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotInst = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotInst = true;
        nFrames = inputFrames(nTimes);
        
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


%% Present Obstruction Maps

disp('Map Presentation');
disp('-----------------');

disp(' ');

if plotInst || plotMean
    
    % Define Plot Limits
    switch obstructData.targetPlane.orientation

        case 'YZ'
            xLimsData = obstructData.targetPlane.position;
            yLimsData = [min(obstructData.positionGrid(:,2)); max(obstructData.positionGrid(:,2))];
            zLimsData = [min(obstructData.positionGrid(:,3)); max(obstructData.positionGrid(:,3))];

        case 'XZ'
            xLimsData = [min(obstructData.positionGrid(:,1)); max(obstructData.positionGrid(:,1))];
            yLimsData = obstructData.targetPlane.position;
            zLimsData = [min(obstructData.positionGrid(:,3)); max(obstructData.positionGrid(:,3))];

        case 'XY'
            xLimsData = [min(obstructData.positionGrid(:,1)); max(obstructData.positionGrid(:,1))];
            yLimsData = [min(obstructData.positionGrid(:,2)); max(obstructData.positionGrid(:,2))];
            zLimsData = obstructData.targetPlane.position;

    end
    
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

    end
    
    if normalise
        
        if contains(caseID, ["Run_Test", "Windsor"])
            xLimsPlot = round((xLimsPlot / 1.044), spacePrecision);
            yLimsPlot = round((yLimsPlot / 1.044), spacePrecision);
            zLimsPlot = round((zLimsPlot / 1.044), spacePrecision);
        end
        
    end
    
    orientation = obstructData.targetPlane.orientation;
    positionData = obstructData.positionGrid;
    mapPerim = [];
    cMap = flipud(viridis(32));
    contourlines = [(0.02 * max(obstructData.mean.areaDensity)); (0.02 * max(obstructData.mean.areaDensity))];
    refPoint = [];
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    nPlanes = 1;
    planeNo = 1;    
end

if plotPerspective
    disp('    Presenting Ray Perspective...');
    
    rayLL = find((obstructData.positionGrid(:,2) == min(obstructData.positionGrid(:,2)) & ...
                 (obstructData.positionGrid(:,3) == min(obstructData.positionGrid(:,3)))));
    rayUL = find((obstructData.positionGrid(:,2) == min(obstructData.positionGrid(:,2)) & ...
                     (obstructData.positionGrid(:,3) == max(obstructData.positionGrid(:,3)))));
    rayLR = find((obstructData.positionGrid(:,2) == max(obstructData.positionGrid(:,2)) & ...
                     (obstructData.positionGrid(:,3) == min(obstructData.positionGrid(:,3)))));
    rayUR = find((obstructData.positionGrid(:,2) == max(obstructData.positionGrid(:,2)) & ...
                     (obstructData.positionGrid(:,3) == max(obstructData.positionGrid(:,3)))));
    rayID = [rayLL; rayUL; rayLR; rayUR]; clear rayLL rayUL rayLR rayUR;
    
    figName = 'Ray_Perspective';
    originPoint = obstructData.originPoint;
    figSubtitle = ' ';
    
    fig = plotVolumeRays(positionData, fig, figName, geometry, originPoint, rayID, figTitle, figSubtitle, ...
                         xLimsPlot, yLimsPlot, zLimsPlot, figSave);
    
    disp(' ');
end

if plotMean
    disp('    Presenting Time-Averaged Obstruction Map...');
    
    scalarData = obstructData.mean.areaDensity;
    figName = 'Time_Averaged_Obstruction_Map';
    figSubtitle = ' ';
    cLims = [0; 0.012]; % cLims = [0; max(scalarData)];
    
    [fig, planeNo] = plotPlanarScalarField(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                           mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, contourlines, ...
                                           xDims, yDims, zDims, refPoint, figTitle, figSubtitle, cLims, ...
                                           xLimsPlot, yLimsPlot, zLimsPlot, normalise, figSave);
    
    disp(' ');
end

if plotInst
    disp('    Presenting Instantaneous Obstruction Map(s)...');
    
    cLims = [0; 0.018]; % cLims = [0; max(cellfun(@max, obstructData.inst.areaDensity))];
    
    figHold = fig;
    
    for i = 1:nFrames
        
        if i ~= 1
            clf(fig);
            fig = figHold;
        end
        
        scalarData = obstructData.inst.areaDensity{i};
        figTime = num2str(obstructData.time(i), ['%.', num2str(timePrecision), 'f']);
        figName = ['Instantaneous_Obstruction_Map_T_', erase(figTime, '.')];
        figSubtitle = [figTime, ' \it{s}'];

        [fig, planeNo] = plotPlanarScalarField(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                               mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, contourlines, ...
                                               xDims, yDims, zDims, refPoint, figTitle, figSubtitle, cLims, ...
                                               xLimsPlot, yLimsPlot, zLimsPlot, normalise, figSave);
    end
    clear i;
    
    disp(' ');
end

if ~plotPerspective && ~plotMean && ~plotInst
    disp('    Skipping Line of Sight Presentation');
    disp(' ');
end

disp(' ');


%% Save Line of Sight Obstruction Data

disp('Data Save Options');
disp('------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Save Data for Future Use? [y/n]: ', 's');
    
    if selection == 'n' | selection == 'N' %#ok<OR2>
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        fileID = erase(['Plane_', obstructData.targetPlane.orientation, '_', ...
                        num2str(obstructData.targetPlane.position), '_Origin_', ...
                        num2str(obstructData.originPoint(1)), '_', ...
                        num2str(obstructData.originPoint(2)), '_', ...
                        num2str(obstructData.originPoint(3))], '.');
        
        obstructData = orderfields(obstructData, {'targetPlane', 'originPoint', 'raysOfInterest', ...
                                                  'positionGrid', 'time', 'inst', 'mean'});
        
        % Save Data
        switch format
            
            case 'A'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/volumeObstruction/', caseID, '/nearField/', dataID], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/volumeObstruction/', caseID, '/nearField/', dataID]);
                end
                
            case 'B'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/volumeObstruction/', caseID, '/farField/', dataID], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/volumeObstruction/', caseID, '/farField/', dataID]);
                end
                
        end
        
        switch format
            
            case 'A'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/volumeObstruction/', caseID, '/nearField/', dataID, '/', fileID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/volumeObstruction/', caseID, '/nearField/', dataID, '/', fileID, '.mat'], ...
                      'caseID', 'dataID', 'fileID', 'obstructData', 'cellSize', 'dL', 'sampleInterval', 'timePrecision', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/volumeObstruction/', caseID, '/farField/', dataID, '/', fileID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/volumeObstruction/', caseID, '/farField/', dataID, '/', fileID, '.mat'], ...
                      'caseID', 'dataID', 'fileID', 'obstructData', 'cellSize', 'dL', 'sampleInterval', 'timePrecision', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp('        Success');
        
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end


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
