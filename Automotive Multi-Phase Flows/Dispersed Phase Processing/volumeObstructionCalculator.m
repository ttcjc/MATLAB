%% Lagrangian Line of Sight Obstruction Calculator v2.2
% ----
% Calculate Mass of Spray Contained Along Rays Projecting From a Point of Interest to an Upstream Plane
% (Based on Volumetric Spray Data Collected Using 'volumeFieldGenerator')


%% Preamble

run preamble;

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
                disp(['Loading ''', fileName, '''...']);
                
                campaignID = load([filePath, fileName], 'campaignID').campaignID;
                caseID = load([filePath, fileName], 'caseID').caseID;
                dataID = load([filePath, fileName], 'dataID').dataID;
                volumeData = load([filePath, fileName], 'volumeData').volumeData;
                cellSize = load([filePath, fileName], 'cellSize').cellSize;
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
            
            if contains(filePath, '/midWake')
                disp(['Loading ''', fileName, '''...']);
                
                campaignID = load([filePath, fileName], 'campaignID').campaignID;
                caseID = load([filePath, fileName], 'caseID').caseID;
                volumeData = load([filePath, fileName], 'volumeData').volumeData;
                cellSize = load([filePath, fileName], 'cellSize').cellSize;
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
            
        case 'C'
            
            if contains(filePath, '/farWake')
                disp(['Loading ''', fileName, '''...']);
                
                campaignID = load([filePath, fileName], 'campaignID').campaignID;
                caseID = load([filePath, fileName], 'caseID').caseID;
                volumeData = load([filePath, fileName], 'volumeData').volumeData;
                cellSize = load([filePath, fileName], 'cellSize').cellSize;
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

disp(' ');
disp(' ');


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry;

disp(' ');
disp(' ');


%% Select Plane of Interest

% Select Plane of Interest
obstructData.obstruction = identifyVolumeSlices(volumeData.positionGrid, spacePrecision, false, normDims);
obstructData.obstruction.planePosition = obstructData.obstruction.position;
obstructData.obstruction = rmfield(obstructData.obstruction, 'position');

% Extract Planar Position Data
index = find(volumeData.positionGrid(:,1) == obstructData.obstruction.planePosition);
obstructData.positionGrid = volumeData.positionGrid(index,:);

nCells = height(obstructData.positionGrid);

disp(' ');
disp(' ');


%% Select Origin Point

disp('Ray Origin Definition');
disp('----------------------');

% Select Ray Origin Point
obstructData.obstruction.pointPosition = zeros([1,3]);

valid = false;
while ~valid
    disp(' ')
    disp('Specify Ray Origin Point:')
    
    obstructData.obstruction.pointPosition(1) = inputPos('X');
    obstructData.obstruction.pointPosition(2) = inputPos('Y');
    obstructData.obstruction.pointPosition(3) = inputPos('Z');

    if (obstructData.obstruction.pointPosition(1) < min(volumeData.positionGrid(:,1)) || obstructData.obstruction.pointPosition(1) > max(volumeData.positionGrid(:,1))) || ...
       (obstructData.obstruction.pointPosition(2) < min(volumeData.positionGrid(:,2)) || obstructData.obstruction.pointPosition(2) > max(volumeData.positionGrid(:,1))) || ...
       (obstructData.obstruction.pointPosition(3) < min(volumeData.positionGrid(:,3)) || obstructData.obstruction.pointPosition(3) > max(volumeData.positionGrid(:,1)))
        disp('        WARNING: Origin Point Lies Outside Volume');
        
        continue;
    end
    
    switch obstructData.obstruction.orientation
        
        case 'YZ'
            
            if obstructData.obstruction.pointPosition(1) < obstructData.obstruction.planePosition
                disp('        WARNING: Origin Point Must Lie in the Positive X-Direction of the Target Plane');
                
                continue;
            end
            
        case 'XZ'
            
            if obstructData.obstruction.pointPosition(2) > obstructData.obstruction.planePosition
                disp('        WARNING: Origin Point Must Lie in the Negative Y-Direction of the Target Plane');
                
                continue;
            end
            
        case 'XY'
            
            if obstructData.obstruction.pointPosition(3) > obstructData.obstruction.planePosition
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

% % Check if Plane of Interest Intersects Geometry
% removeIntersect = false;
% 
% switch obstructData.obstruction.orientation
%     
%     case 'YZ'
%         
%         if obstructData.obstruction.planePosition <= xDims(2)
%             removeIntersect = true;
%         end
%         
%     case 'XZ'
%         
%         if obstructData.obstruction.planePosition >= yDims(1)
%             removeIntersect = true;
%         end
%         
%     case 'XY'
%         
%         if obstructData.obstruction.planePosition <= zDims(2)
%             removeIntersect = true;
%         end
%         
% end
% 
% % Remove Erroneous Data From Cells Intersecting Geometry
% if removeIntersect
%     disp('        Removing Erroneous Data From Grid Cells Intersecting Geometry...');
% 
%     % Perform Removal
%     volumeDataVars = fieldnames(volumeData);
%     nonFieldVars = {'positionGrid'; 'time'};
%     fieldVars = setdiff(volumeDataVars, nonFieldVars);
%     clear volumeDataVars nonFieldVars;
%     
%     parts = fieldnames(geometry);
%     for i = 1:height(parts)
%         geoPoints = unique(geometry.(parts{i}).vertices, 'rows');
%         DT = delaunay(geoPoints);
% 
%         index = ~isnan(tsearchn(geoPoints, DT, volumeData.positionGrid));
% 
%         for j = 1:height(fields)
%             volumeData.(fieldVars{j}).mean(index,:) = NaN;
% 
%             for k = 1:nTimes
%                 volumeData.(fieldVars{j}).inst{k}(index,:) = NaN;
%             end
% 
%         end
%         clear j;
% 
%     end
%     clear i parts fieldVars;
%     
% end

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

% Small Distance Used During LOS Numerical Integration
dL = 1.25e-4;

if ~normDims
    
    if strcmp(campaignID, 'Windsor_fullScale')
        dL = 5e-4;
    elseif strcmp(campaignID, 'Windsor_Upstream_2023')
        dL = 1.25e-4;
    end
    
end

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Obstruction', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);
parforWaitBar(wB, nTimes);

% Perform Calculation
areaDensity = cell(nTimes,1); areaDensity(:) = {zeros([nCells,1])};

density = volumeData.density.inst;
positionGrid = obstructData.positionGrid;
pointPosition = obstructData.obstruction.pointPosition;
parfor i = 1:nTimes
    densityField = reshape(full(density{i}), gridShape);
    densityInterp = griddedInterpolant(x, y, z, densityField, 'linear', 'none');
    
    for j = 1:height(positionGrid)
        dirVec = positionGrid(j,:) - pointPosition;
        distFull = sqrt(dirVec(1)^2 + dirVec(2)^2 + dirVec(3)^2);

        dist = (dL:dL:distFull)';
        samplePoints = pointPosition + (dist * (dirVec / distFull));
        
        areaDensity{i}(j) = sum((densityInterp(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3))) * dL);
    end
    
    % Make Arrays Sparse
    areaDensity{i} = sparse(areaDensity{i});
    
    % Remove Unnecessary Data
    density{i} = [];
    
    % Update Waitbar
    send(dQ, []);
end
clear density positionGrid pointPosition;

delete(wB);

obstructData.areaDensity.inst = areaDensity; clear areaDensity;

disp(' ');

disp('    Calculating Time-Averaged Obstruction...');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Time-Averaged Obstruction', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Perform Calculation
obstructData.areaDensity.mean = zeros([nCells,1]);

densityField = reshape(full(volumeData.density.mean), gridShape);
densityInterp = griddedInterpolant(x, y, z, densityField, 'linear', 'none');

for i = 1:height(obstructData.positionGrid)
    dirVec = obstructData.positionGrid(i,:) - obstructData.obstruction.pointPosition;
    distFull = sqrt(dirVec(1)^2 + dirVec(2)^2 + dirVec(3)^2);

    dist = (dL:dL:distFull)';
    samplePoints = obstructData.obstruction.pointPosition + (dist * (dirVec / distFull));
    
    obstructData.areaDensity.mean(i) = sum((densityInterp(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3))) * dL);
    
    % Update Waitbar
    waitbar((i / height(obstructData.positionGrid)), wB);
end
clear i;

obstructData.areaDensity.mean = sparse(obstructData.areaDensity.mean);

delete(wB);

clear volumeData;

obstructID = erase(['Plane_', obstructData.obstruction.orientation, '_', ...
                num2str(obstructData.obstruction.planePosition), '_Origin_', ...
                num2str(obstructData.obstruction.pointPosition(1)), '_', ...
                num2str(obstructData.obstruction.pointPosition(2)), '_', ...
                num2str(obstructData.obstruction.pointPosition(3))], '.');

obstructData = orderfields(obstructData, {'positionGrid', 'time', 'areaDensity', ...
                                          'obstruction'});

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


%% Present Obstruction Maps

disp('Map Presentation');
disp('-----------------');

disp(' ');

if plotMean || plotInst
    
    parts = fieldnames(geometry);
    for i = 1:height(parts)
        geometry.(parts{i}).vertices = geometry.(parts{i}).vertices / normLength;
    end
    clear i parts;
    
    orientation = obstructData.obstruction.orientation;
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
    
    switch orientation

        case 'YZ'
            xLimsData = obstructData.obstruction.planePosition;
            yLimsData = [min(obstructData.positionGrid(:,2)); max(obstructData.positionGrid(:,2))];
            zLimsData = [min(obstructData.positionGrid(:,3)); max(obstructData.positionGrid(:,3))];

        case 'XZ'
            xLimsData = [min(obstructData.positionGrid(:,1)); max(obstructData.positionGrid(:,1))];
            yLimsData = obstructData.obstruction.planePosition;
            zLimsData = [min(obstructData.positionGrid(:,3)); max(obstructData.positionGrid(:,3))];

        case 'XY'
            xLimsData = [min(obstructData.positionGrid(:,1)); max(obstructData.positionGrid(:,1))];
            yLimsData = [min(obstructData.positionGrid(:,2)); max(obstructData.positionGrid(:,2))];
            zLimsData = obstructData.obstruction.planePosition;

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

if plotMean
    disp('    Presenting Time-Averaged Obstruction Map...');
    
    scalarData = full(obstructData.areaDensity.mean);
    figName = ['Average_Obstruction_', caseID, '_', obstructID];
    
    if strcmp(campaignID, 'Windsor_fullScale')
        
        if contains(caseID, 'uncoupled')
            contourlines = [0.02; 0.02] * 2.45;
        else
%             contourlines = [0.02; 0.02] * 6.25; % 1.75L
            contourlines = [0.02; 0.02] * 11.5; % 3.75L
        end
        
    elseif strcmp(campaignID, 'Windsor_Upstream_2023')
        contourlines = [0.02; 0.02] * 7.2e-3;
    end
    
    if ~exist('contourlines', 'var')
        contourlines = [];
    end
    
    figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
    
    if strcmp(campaignID, 'Windsor_fullScale')

        if contains(caseID, 'uncoupled')
            cLims = [0; 2.45];
        else
%             cLims = [0; 6.25];  % 1.75L
            cLims = [0; 11.5];  % 3.75L
        end

    elseif strcmp(campaignID, 'Windsor_Upstream_2023')
        cLims = [0; 7.2e-3];
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
    disp('    Presenting Instantaneous Obstruction Map(s)...');
    
    contourlines = [];
    
    if strcmp(campaignID, 'Windsor_fullScale')

        if contains(caseID, 'uncoupled')
            cLims = [0; 4.25];
        else
%             cLims = [0; 7.6]; % 1.75
            cLims = [0; 13.65]; % 3.75
        end

    elseif strcmp(campaignID, 'Windsor_Upstream_2023')
        cLims = [0; 14e-3];
    end
    
    if ~exist('cLims', 'var')
        instMax = 0;

        for j = 1:nTimes
            instMax = max(instMax, max(obstructData.areaDensity.inst{j}));
        end

        cLims = full([0; instMax]);
    end
    
    figHold = fig;
    
    for i = startFrame:endFrame
        
        if i ~= startFrame
            clf(fig);
            fig = figHold;
        end
        
        scalarData = full(obstructData.areaDensity.inst{i});
        figTime = num2str(obstructData.time(i), ['%.', num2str(timePrecision), 'f']);
        figName = ['Inst_Obstruction_T', erase(figTime, '.'), '_', caseID, '_', obstructID];
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

if ~plotMean && ~plotInst
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
        
        % Save Data
        switch format
            
            case 'A'
                
                if ~exist([saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/nearWake/', dataID], 'dir')
                    mkdir([saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/nearWake/', dataID]);
                end
                
            case 'B'
                
                if ~exist([saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/midWake/', dataID], 'dir')
                    mkdir([saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/midWake/', dataID]);
                end
                
            case 'C'
                
                if ~exist([saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/farWake/', dataID], 'dir')
                    mkdir([saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/farWake/', dataID]);
                end
                
        end
        
        switch format
            
            case 'A'
                disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/nearWake/', dataID, '/', obstructID, '.mat']);
                
                save([saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/nearWake/', dataID, '/', obstructID, '.mat'], ...
                     'campaignID', 'caseID', 'dataID', 'obstructID', 'obstructData', 'cellSize', 'dL', 'sampleInt', 'timePrecision', 'dLims', 'normDims', '-v7.3', '-noCompression');
                
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/midWake/', dataID, '/', obstructID, '.mat']);
                
                save([saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/midWake/', dataID, '/', obstructID, '.mat'], ...
                      'campaignID',  'caseID', 'dataID', 'obstructID', 'obstructData', 'cellSize', 'dL', 'sampleInt', 'timePrecision', 'dLims', 'normDims', '-v7.3', '-noCompression');
                
                disp('        Success');
                
            case 'C'
                disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/farWake/', dataID, '/', obstructID, '.mat']);
                
                save([saveLoc, '/Numerical/MATLAB/volumeObstruction/', campaignID, '/', caseID, '/farWake/', dataID, '/', obstructID, '.mat'], ...
                      'campaignID',  'caseID', 'dataID', 'obstructID', 'obstructData', 'cellSize', 'dL', 'sampleInt', 'timePrecision', 'dLims', 'normDims', '-v7.3', '-noCompression');
                
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


function frameNo = inputFrames(Nt, type)

    frameNo = str2double(input(['    Input Desired ', type, ' Frame [1-', num2str(Nt), ']: '], 's'));
    
    if isnan(frameNo) || frameNo < 1 || frameNo > Nt
        disp('        WARNING: Invalid Entry');
        
        frameNo = -1;
    end

end
