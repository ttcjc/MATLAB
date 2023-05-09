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

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallelisation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures


%% Lagrangian Line of Sight Calculator v1.0

dL = 1.25e-4; % Small Distance Used During LoS Numerical Integration

% samplesPerCell = 32; % Number of Samples per Cell Used During LoS Interpolation

disp('=============================');
disp('Line of Sight Calculator v1.1');
disp('=============================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - 


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
LoSdata = identifyVolumeSlices(volumeData.positionGrid, spacePrecision, false);

% Extract Planar Position Data
index = find(volumeData.positionGrid(:,1) == LoSdata.planeLocation);
LoSdata.positionGrid = volumeData.positionGrid(index,:);

disp(' ');
disp(' ');


%% Select Origin Point

disp('Origin Point Definition');
disp('------------------------');

% Select Origin Point
LoSdata.originPoint = zeros(1,3);

valid = false;
while ~valid
    disp(' ')
    disp('Specify Origin Point:')
    
    LoSdata.originPoint(1) = inputPos('X');
    LoSdata.originPoint(2) = inputPos('Y');
    LoSdata.originPoint(3) = inputPos('Z');

    if (LoSdata.originPoint(1) < min(volumeData.positionGrid(:,1)) || LoSdata.originPoint(1) > max(volumeData.positionGrid(:,1))) || ...
       (LoSdata.originPoint(2) < min(volumeData.positionGrid(:,2)) || LoSdata.originPoint(2) > max(volumeData.positionGrid(:,1))) || ...
       (LoSdata.originPoint(3) < min(volumeData.positionGrid(:,3)) || LoSdata.originPoint(3) > max(volumeData.positionGrid(:,1)))
        disp('        WARNING: Origin Point Lies Outside Volume');
        continue;
    end
    
    switch LoSdata.orientation
        
        case 'YZ'
            
            if LoSdata.originPoint(1) < LoSdata.planeLocation
                disp('        WARNING: Origin Point Must Lie in the Positive X-Direction of the Target Plane');
                continue;
            end
            
        case 'XZ'
            
            if LoSdata.originPoint(2) > LoSdata.planeLocation
                disp('        WARNING: Origin Point Must Lie in the Negative Y-Direction of the Target Plane');
                continue;
            end
            
        case 'XY'
            
            if LoSdata.originPoint(3) > LoSdata.planeLocation
                disp('        WARNING: Origin Point Must Lie in the Positive Z-Direction of the Target Plane');
                continue;
            end
            
    end
    
    valid = true;    
end

% Restore Normalised Dimensions
if normalise

    if contains(caseName, ["Run_Test", "Windsor"])
        volumeData.positionGrid = round((volumeData.positionGrid / 1.044), spacePrecision);
        LoSdata.originPoint = round((LoSdata.originPoint / 1.044), spacePrecision);
        LoSdata.planeLocation = round((LoSdata.planeLocation / 1.044), spacePrecision);
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

switch LoSdata.orientation
    
    case 'YZ'
        
        if LoSdata.planeLocation <= xDims(2)
            removeIntersect = true;
        end
        
    case 'XZ'
        
        if LoSdata.planeLocation >= yDims(1)
            removeIntersect = true;
        end
        
    case 'XY'
        
        if LoSdata.planeLocation <= zDims(2)
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

LoSdata.inst.time = volumeData.inst.time;

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Line of Sight', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, height(LoSdata.inst.time));

% Calculate Instantaneous Line of Sight
massLoS = cell(height(LoSdata.inst.time),1);
nParticlesLoS = massLoS;
dLoS = massLoS;
mass = volumeData.inst.mass;
nParticles = volumeData.inst.nParticles;
d = volumeData.inst.d10;
cellVolume = cellSize.volume;
positionGrid = LoSdata.positionGrid;
% orientation = LoSdata.orientation;
% cellSizeX = cellSize.x;
% cellSizeY = cellSize.y;
% cellSizeZ = cellSize.z;
cellSizeTarget = cellSize.target;
% originPointX = LoSdata.originPoint(1);
% originPointY = LoSdata.originPoint(2);
% originPointZ = LoSdata.originPoint(3);
originPoint = LoSdata.originPoint;
% planeLocation = LoSdata.planeLocation;
parfor i = 1:height(LoSdata.inst.time)
    massDensity = reshape((mass{i} / cellVolume), gridShape);
    nParticleDensity = reshape((nParticles{i} / cellVolume), gridShape);
    dDensity = reshape((d{i} / cellVolume), gridShape);

    massInterp = griddedInterpolant(x, y, z, massDensity, 'linear', 'none');
    nParticlesInterp = griddedInterpolant(x, y, z, nParticlesDensity, 'linear', 'none');
    dInterp = griddedInterpolant(x, y, z, dDensity, 'linear', 'none');
    
    massLoS{i} = zeros(height(positionGrid),1);
    nParticlesLoS{i} = massLoS{i};
    dLoS{i} = massLoS{i};
    
%     switch orientation
%         
%         case 'YZ'
%             
%             dX = cellSizeX / samplesPerCell;
%             sampleX = (originPointX:-dX:planeLocation)';
%             for j = 1:height(positionGrid)
%                 t = (sampleX - originPointX) / (positionGrid(j,1) - originPointX);
%                 sampleXYZ = ([originPointX, originPointY, originPointZ] + t*(positionGrid(j,:) - [originPointX, originPointY, originPointZ]));
%                 interpMass = interp(sampleXYZ(:,1), sampleXYZ(:,2), sampleXYZ(:,3)) * dX;
%                 
%                 if sum(isnan(interpMass)) > samplesPerCell
%                     massInLoS{i}(j) = NaN;
%                 else
%                     massInLoS{i}(j) = sum(interpMass, 'omitNaN');
%                 end
%                 
%             end
%         
%         case 'XZ'
% 
%             dY = cellSizeY / samplesPerCell;
%             sampleY = (originPointY:dY:planeLocation)';
%             for j = 1:height(positionGrid)
%                 t = (sampleY -originPointY) / (positionGrid(j,2) - originPointY);
%                 sampleXYZ = ([originPointX, originPointY, originPointZ] + t*(positionGrid(j,:) - [originPointX, originPointY, originPointZ]));
%                 interpMass = interp(sampleXYZ(:,1), sampleXYZ(:,2), sampleXYZ(:,3)) * dY;
%                 
%                 if sum(isnan(interpMass)) > samplesPerCell
%                     massInLoS{i}(j) = NaN;
%                 else
%                     massInLoS{i}(j) = sum(interpMass, 'omitNaN');
%                 end
%                 
%             end
%             
%         case 'XY'
% 
%             dZ = cellSizeZ / samplesPerCell;
%             sampleZ = (originPointZ:-dZ:planeLocation)';
%             for j = 1:height(positionGrid)
%                 t = (sampleZ - originPointZ) / (positionGrid(j,3) - originPointZ);
%                 sampleXYZ = ([originPointX, originPointY, originPointZ] + t*(positionGrid(j,:) - [originPointX, originPointY, originPointZ]));
%                 interpMass = interp(sampleXYZ(:,1), sampleXYZ(:,2), sampleXYZ(:,3)) * dZ;
%                 
%                 if sum(isnan(interpMass)) > samplesPerCell
%                     massInLoS{i}(j) = NaN;
%                 else
%                     massInLoS{i}(j) = sum(interpMass, 'omitNaN');
%                 end
%                 
%             end
%     end
    
    for j = 1:height(positionGrid)
        dirVec = originPoint - positionGrid(j,:);
        distFull = sqrt(dirVec(1)^2 + dirVec(2)^2 + dirVec(3)^2);

        dist = (0:dL:distFull)';
        samplePoints = positionGrid(j,:) + (dist * (dirVec / distFull));

        massOnLine = massInterp(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3)) * dL;
        nParticlesOnLine = nParticlesInterp(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3)) * dL;
        dOnLine = dInterp(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3)) * dL;
        
        if sum(isnan(massOnLine)) > ceil(2 * (cellSizeTarget / dL))
            massLoS{i}(j) = NaN;
            nParticleLoS{i}(j) = NaN;
            dLoS{i}(j) = NaN;
        else
            massLoS{i}(j) = sum(massOnLine, 'omitNaN');
            nParticlesLoS{i}(j) = sum(nParticlesOnLine, 'omitNaN');
            dLoS{i}(j) = sum(dOnLine, 'omitNaN');
        end

    end
    
    send(dQ, []);
end
LoSdata.inst.mass = massLoS;
% clear massLoS mass cellVolume positionGrid orientation cellSizeX cellSizeY cellSizeZ originPointX originPointY originPointZ planeLocation;
clear massLoS nParticlesLoS dLoS mass nParticles d cellVolume positionGrid cellSizeTarget originPoint;

delete(wB);

disp(' ');

disp('    Calculating Time-Averaged Line of Sight...');

% Calculate Time-Averaged Line of Sight
massDensity = reshape((volumeData.mean.mass / cellSize.volume), gridShape);
massInterp = griddedInterpolant(x, y, z, massDensity, 'linear', 'none');

LoSdata.mean.mass = zeros(height(LoSdata.positionGrid),1);

% switch LoSdata.orientation
% 
%     case 'YZ'
% 
%         dX = cellSize.x / samplesPerCell;
%         sampleX = (LoSdata.originPoint(1):-dX:LoSdata.planeLocation)';
%         for i = 1:height(LoSdata.positionGrid)
%             t = (sampleX - LoSdata.originPoint(1)) / (LoSdata.positionGrid(i,1) - LoSdata.originPoint(1));
%             sampleXYZ = (LoSdata.originPoint + t*(LoSdata.positionGrid(i,:) - LoSdata.originPoint));
%             interpMass = interp(sampleXYZ(:,1), sampleXYZ(:,2), sampleXYZ(:,3)) * dX;
%             
%             if sum(isnan(interpMass)) > samplesPerCell
%                 LoSdata.mean.mass(i) = NaN;
%             else
%                 LoSdata.mean.mass(i) = sum(interpMass, 'omitNaN');
%             end
%             
%         end
% 
%     case 'XZ'
% 
%         dY = cellSize.y / samplesPerCell;
%         sampleY = (LoSdata.originPoint(2):dY:LoSdata.planeLocation)';
%         for i = 1:height(LoSdata.positionGrid)
%             t = (sampleY - LoSdata.originPoint(2)) / (LoSdata.positionGrid(i,2) - LoSdata.originPoint(2));
%             sampleXYZ = (LoSdata.originPoint + t*(LoSdata.positionGrid(i,:) - LoSdata.originPoint));
%             interpMass = interp(sampleXYZ(:,1), sampleXYZ(:,2), sampleXYZ(:,3)) * dY;
%             
%             if sum(isnan(interpMass)) > samplesPerCell
%                 LoSdata.mean.mass(i) = NaN;
%             else
%                 LoSdata.mean.mass(i) = sum(interpMass, 'omitNaN');
%             end
%             
%         end
% 
%     case 'XY'
% 
%         dZ = cellSize.z / samplesPerCell;
%         sampleZ = (LoSdata.originPoint(3):-dZ:LoSdata.planeLocation)';
%         for i = 1:height(LoSdata.positionGrid)
%             t = (sampleZ - LoSdata.originPoint(3)) / (LoSdata.positionGrid(i,3) - LoSdata.originPoint(3));
%             sampleXYZ = (LoSdata.originPoint + t*(LoSdata.positionGrid(i,:) - LoSdata.originPoint));
%             interpMass = interp(sampleXYZ(:,1), sampleXYZ(:,2), sampleXYZ(:,3)) * dZ;
%             
%             if sum(isnan(interpMass)) > samplesPerCell
%                 LoSdata.mean.mass(i) = NaN;
%             else
%                 LoSdata.mean.mass(i) = sum(interpMass, 'omitNaN');
%             end
%         
%         end
% 
% end

for i = 1:height(LoSdata.positionGrid)
    dirVec = LoSdata.originPoint - LoSdata.positionGrid(i,:);
    distFull = sqrt(dirVec(1)^2 + dirVec(2)^2 + dirVec(3)^2);

    dist = (0:dL:distFull)';
    samplePoints = LoSdata.positionGrid(i,:) + (dist * (dirVec / distFull));

    massOnLine = massInterp(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3)) * dL;
    
    if sum(isnan(massOnLine)) > ceil(2 * (cellSize.target / dL))
        LoSdata.mean.mass(i) = NaN;
    else
        LoSdata.mean.mass(i) = sum(massOnLine, 'omitNaN');
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

valid = false;
while ~valid
    disp(' ');
    selection = input('Plot Instantaneous Maps? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotInst = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotInst = true;
        nFrames = inputFrames(height(LoSdata.inst.time));
        
        if nFrames == -1
            continue;
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end

disp(' ');
disp(' ');


%% Present Contaminant Maps

disp('Map Presentation');
disp('-----------------');

disp(' ');

if plotInst || plotMean
    
    switch LoSdata.orientation

        case 'YZ'
            xLimsData = LoSdata.planeLocation;
            yLimsData = [min(LoSdata.positionGrid(:,2)); max(LoSdata.positionGrid(:,2))];
            zLimsData = [min(LoSdata.positionGrid(:,3)); max(LoSdata.positionGrid(:,3))];

        case 'XZ'
            xLimsData = [min(LoSdata.positionGrid(:,1)); max(LoSdata.positionGrid(:,1))];
            yLimsData = LoSdata.planeLocation;
            zLimsData = [min(LoSdata.positionGrid(:,3)); max(LoSdata.positionGrid(:,3))];

        case 'XY'
            xLimsData = [min(LoSdata.positionGrid(:,1)); max(LoSdata.positionGrid(:,1))];
            yLimsData = [min(LoSdata.positionGrid(:,2)); max(LoSdata.positionGrid(:,2))];
            zLimsData = LoSdata.planeLocation;

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
            end

    end

    if normalise
        xLimsPlot = round((xLimsPlot / 1.044), spacePrecision);
        yLimsPlot = round((yLimsPlot / 1.044), spacePrecision);
        zLimsPlot = round((zLimsPlot / 1.044), spacePrecision);
    end
    
    orientation = LoSdata.orientation;
    positionData = LoSdata.positionGrid;
    mapPerim = [];
    cMap = flipud(viridis(32));
    contourlines = [];
    CoM = LoSdata.originPoint;
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    nPlanes = 1;
    planeNo = 1;
    
end

if plotMean
    disp('    Presenting Time-Averaged Line of Sight Map...');
    
    scalarData = LoSdata.mean.mass;
    figName = 'Time_Averaged_LoS_Map';
    figSubtitle = ' ';
    cLims = [0; max(LoSdata.mean.mass)];
%     cLims = [0; 5e-3];
    
    [fig, planeNo] = plotPlanarScalarField(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                           mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, contourlines, ...
                                           xDims, yDims, zDims, CoM, figTitle, figSubtitle, cLims, ...
                                           xLimsPlot, yLimsPlot, zLimsPlot, normalise);
    
    disp(' ');
end

if plotInst
    disp('    Presenting Instantaneous Line of Sight Map(s)...');
    
    cLims = [0; max(cellfun(@max, LoSdata.inst.mass))];
%     cLims = [0; 5e-3];
    
    figHold = fig;
    
    for i = 1:nFrames
        
        if i ~= 1
            clf(fig);
            fig = figHold;
        end
        
        scalarData = LoSdata.inst.mass{i};
        figTime = num2str(LoSdata.inst.time(i), ['%.', num2str(timePrecision), 'f']);
        figName = ['Instantaneous_LoS_Map_T_', erase(figTime, '.')];
        figSubtitle = [figTime, ' \it{s}'];

        [fig, planeNo] = plotPlanarScalarField(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                               mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, contourlines, ...
                                               xDims, yDims, zDims, CoM, figTitle, figSubtitle, cLims, ...
                                               xLimsPlot, yLimsPlot, zLimsPlot, normalise);
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