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

nProc = 4; % Number of Processors Used for Process-Based Parallelisation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures


%% Lagrangian Volume Field Generator v4.0

cloudName = 'kinematicCloud'; % OpenFOAM Cloud Name

figSave = false; % Save .fig File(s)

normalise = false; % Normalisation of Dimensions

disp('===========================');
disp('Volume Field Generator v4.0');
disp('===========================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Updated calls to 'globalPos' to 'positionCartesian'
% v1.2 - Updated to Support Changes to 'timeDirectories.m'
% v2.0 - Rewritten to Follow Recent Lagrangian Processing Structure Changes
% v2.1 - Added Time-Averaging Functionality
% v3.0 - Rewrite, Accommodating New OpenFOAM Data Formats
% v3.1 - Added Support for Arithmetic and Sauter Mean Particle Diameters
% v3.2 - Changed to Thread-Based Parallelization to Reduce Memory Requirements
% v3.3 - Added Support for Full-Scale Windsor Model Simulations
% v3.4 - Changed Primary Output From Mass to Density
% v4.0 - Rewrite, Making Use of Sparse Arrays to Reduce Memory Requirements


%% Initialise Case

[caseFolder, caseID, timeDirs, deltaT, timePrecision, geometry, ...
 xDims, yDims, zDims, spacePrecision, normalise, normLength] = initialiseCaseData(normalise);

disp(' ');
disp(' ');


%% Select Region of Interest

disp('Region of Interest');
disp('-------------------');

disp(' ');

disp('Possible Regions of Interest:');
disp('    A: Near Wake');
disp('    B: Mid Wake');
disp('    C: Far Wake');

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


%% Initialise Lagrangian Data

maxNumCompThreads(nProc);

[dataID, LagProps, ~, ~, LagData, sampleInterval] = initialiseLagData(saveLocation, caseFolder, caseID, ...
                                                                      cloudName, false, false, ...
                                                                      true, timeDirs, deltaT, ...
                                                                      timePrecision, maxNumCompThreads);

if ~contains(caseID, 'Windsor_SB_fullScale_multiPhase')
    disp(' ');
    disp('WARNING: Far-Wake Data Is Unavailable for This Case');
    disp('         Performing Analysis on Mid-Wake Data Instead');
    
    format = 'B';
end
    

disp(' ');
disp(' ');


%% Select Field Options

disp('Field Options');
disp('--------------');

if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
    dLimsDefault = [1; 147];
else
    dLimsDefault = [20; 400];
end

valid = false;
while ~valid
    disp(' ');
    selection = input('Filter Particle Diameters? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        dLims = dLimsDefault;

        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        dLims = zeros([2,1]);
        
        dLims(1) = inputD('Min');

        if dLims(1) == -1
            continue;
        end

        dLims(2) = inputD('Max');

        if dLims(2) == -1
            continue;
        end
        
        dLims = sort(dLims);
        dLims(1) = floor(dLims(1));
        dLims(2) = ceil(dLims(2));
        
        if dLims(2) < dLimsDefault(1) || dLims(1) > dLimsDefault(2)
            disp('        WARNING: No Lagrangian Data in Diameter Range');
            continue;
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid dLimsDefault;

% Generate Dataset ID
if normalise
    dataID = [dataID, '_D', num2str(dLims(1)), '_D', num2str(dLims(2)), '_Norm'];
else
    dataID = [dataID, '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];
end

disp(' ');
disp(' ');


%% Generate Volume Field

disp('Volume Field Generation');
disp('------------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

evalc('parpool(''threads'');');

%%%%

disp(' ');

disp('    Initialising...');

nTimes = height(LagData.time);

% Identify Empty Time Instances Not Related to Case Initialisation
emptyTimes = cellfun(@isempty, LagData.origId);
firstValidTime = find(emptyTimes == false, 1, 'first');

suspectTimes = find(emptyTimes(firstValidTime:end) == true) + (firstValidTime - 1);

if ~isempty(suspectTimes)
    
    for i = 1:height(suspectTimes)
        disp(['        WARNING: Time ''', num2str(LagData.time(suspectTimes(i))), ''' Is Unexpectedly Empty']);
    end
    clear i;
    
end

% Shift Data Origin
if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
    
    for i = 1:nTimes
        
        if ~isempty(LagData.positionCartesian{i})
            LagData.positionCartesian{i}(:,1) = LagData.positionCartesian{i}(:,1) + 1.325;
        end
        
    end
    clear i;
    
end

% Normalise Coordinate System
if normalise
    disp('        Normalising Coordinate System');
    
    % Initialise Progress Bar
    wB = waitbar(0, 'Normalising Coordinate System', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';

    % Perform Normalisation
    for i = 1:nTimes

        if ~isempty(LagData.positionCartesian{i})
            LagData.positionCartesian{i}  = round((LagData.positionCartesian{i} / normLength), ...
                                                  spacePrecision);
        end

        % Update Waitbar
        waitbar((i / nTimes), wB);
    end
    clear i;

    delete(wB);
end

% Specify Region Boundaries
switch format
    
    case 'A' % 0.75 L
        
        if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
            xLimsData = [0.31875; 1.26625];
            yLimsData = [-0.4176; 0.4176];
            zLimsData = [0; 0.4176];
        else
            xLimsData = [1.275; 5.065];
            yLimsData = [-1.6704; 1.6704];
            zLimsData = [0; 1.6704];
        end
        
    case 'B' % 2 L
        
        if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
            xLimsData = [0.31875; 2.57125];
            yLimsData = [-0.522; 0.522];
            zLimsData = [0; 0.522];
        else
            xLimsData = [1.275; 10.285];
            yLimsData = [-2.088; 2.088];
            zLimsData = [0; 2.088];
        end
        
    case 'C' % 4 L
        
        xLimsData = [1.275; 18.637];
        yLimsData = [-2.5056; 2.5056];
        zLimsData = [0; 2.5056];

end

if normalise
    xLimsData = round((xLimsData / normLength), spacePrecision);
    yLimsData = round((yLimsData / normLength), spacePrecision);
    zLimsData = round((zLimsData / normLength), spacePrecision);
end

disp(' ');

% Collate Particles of Interest
disp('    Collating Particles of Interest...');

% Initialise Progress Bar
wB = waitbar(0, 'Collating Particles of Interest', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);
parforWaitBar(wB, nTimes);

% Perform Collation
index = cell(nTimes,1);

d = LagData.d;
positionCartesian = LagData.positionCartesian;
parfor i = 1:nTimes
    
    if ~isempty(positionCartesian{i})
        index{i} = find(((d{i} * 1e6) >= dLims(1)) & ...
                        ((d{i} * 1e6) <= dLims(2)) & ...
                        (positionCartesian{i}(:,1) >= xLimsData(1)) & ...
                        (positionCartesian{i}(:,1) <= xLimsData(2)) & ...
                        (positionCartesian{i}(:,2) >= yLimsData(1)) & ...
                        (positionCartesian{i}(:,2) <= yLimsData(2)) & ...
                        (positionCartesian{i}(:,3) >= zLimsData(1)) & ...
                        (positionCartesian{i}(:,3) <= zLimsData(2))); %#ok<PFBNS>
    end
    
    % Remove Unnecessary Data
    d{i} = [];
    positionCartesian{i} = [];
    
    % Update Waitbar
    send(dQ, []);
end
clear d positionCartesian;
        
delete(wB);

% Remove Unnecessary Data
disp('        Removing Unnecessary Data...');

LagFields = fieldnames(LagData);
reqFields = {'time'; 'd'; 'nParticle'; 'positionCartesian'};

LagData = rmfield(LagData, LagFields(~ismember(LagFields, reqFields)));
LagProps = LagProps(ismember(LagProps, reqFields));

for i = 1:nTimes
    
    for j = 1:height(LagProps)
        LagData.(LagProps{j}){i} = LagData.(LagProps{j}){i}(index{i},:);
    end
    clear j;

end
clear i;

disp(' ');

% Generate Instantaneous Volume Field
disp('    Generating Presentation Grid...');

% Set Target Spatial Resolution
if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
    cellSize.target = 8e-3;
else
    cellSize.target = 32e-3;
end

% Adjust Uniform Cell Size to Fit Region of Interest
if normalise
    cellSize.target = round((cellSize.target / normLength), spacePrecision);
end

cellSize.x = (xLimsData(2) - xLimsData(1)) / round(((xLimsData(2) - xLimsData(1)) / cellSize.target));
cellSize.y = (yLimsData(2) - yLimsData(1)) / round(((yLimsData(2) - yLimsData(1)) / cellSize.target));
cellSize.z = (zLimsData(2) - zLimsData(1)) / round(((zLimsData(2) - zLimsData(1)) / cellSize.target));

cellSize.volume = cellSize.x * cellSize.y * cellSize.z;

[x, y, z] = ndgrid(xLimsData(1):(cellSize.x):xLimsData(2), ...
                   yLimsData(1):(cellSize.y):yLimsData(2), ...
                   zLimsData(1):(cellSize.z):zLimsData(2));

volumeData.positionGrid = [x(:), y(:), z(:)]; clear x y z;

nCells = height(volumeData.positionGrid);

% Assign Particles to Mesh Nodes
disp('        Assigning Particles to Grid Cells');

% Initialise Progress Bar
wB = waitbar(0, 'Assigning Particles to Grid Cells', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);
parforWaitBar(wB, nTimes);

% Perform Assignment
totalParticles = cellfun(@height, LagData.positionCartesian);
index = cell(nTimes,1); % Array Position of Nearest Mesh Node

xVals = unique(volumeData.positionGrid(:,1));
yVals = unique(volumeData.positionGrid(:,2));
zVals = unique(volumeData.positionGrid(:,3));
positionCartesian = LagData.positionCartesian;
positionGrid = volumeData.positionGrid;
parfor i = 1:nTimes
    
    if totalParticles(i) > 0
        index3D = zeros([totalParticles(i),3], 'uint32');
        positions = zeros([totalParticles(i),3]);
        
        for j = 1:totalParticles(i)
            [~, index3D(j,1)] = min(abs(positionCartesian{i}(j,1) - xVals));
            [~, index3D(j,2)] = min(abs(positionCartesian{i}(j,2) - yVals));
            [~, index3D(j,3)] = min(abs(positionCartesian{i}(j,3) - zVals));
            
            positions(j,:) = [xVals(index3D(j,1)), yVals(index3D(j,2)), zVals(index3D(j,3))];
        end
        
        [~, index{i}] = ismember(positions, positionGrid, 'rows');
    end
    
    % Remove Unnecessary Data
    positionCartesian{i} = [];
    
    % Update Waitbar
    send(dQ, []);
end
clear xVals yVals zVals positionCartesian positionGrid;

delete(wB);

disp(' ');

% Generate Instantaneous Volume Field
disp('    Generating Instantaneous Volume Field...');

volumeData.time = LagData.time;

% Calculate Instantaneous Field Variables
disp('        Calculating Instantaneous Field Variables');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Field Variables', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);
parforWaitBar(wB, nTimes);

% Perform Calculation
nParticles = cell(nTimes,1); nParticles(:) = {zeros([nCells,1])}; % Number of Particles in Cell
density = nParticles; % Spray Density in Cell
d32 = nParticles; % Sauter Mean Diameter in Cell
d10 = nParticles; % Arithmetic Mean Diameter in Cell

d_tmp = zeros([nCells,1]);
nParticle = cellfun(@double, LagData.nParticle, 'uniformOutput', false);
d = cellfun(@double, LagData.d, 'uniformOutput', false);
cellVolume = cellSize.volume;
clear LagData;
parfor i = 1:nTimes
    
    if totalParticles(i) > 0
        d30 = d_tmp;
        d20 = d_tmp;
        
        for j = 1:totalParticles(i)
            nParticles{i}(index{i}(j)) = nParticles{i}(index{i}(j)) + ...
                                         nParticle{i}(j);
            
            density{i}(index{i}(j)) = density{i}(index{i}(j)) + ...
                                      (nParticle{i}(j) * ((1 / 12) * tau * (d{i}(j)^3)));
            
            d30(index{i}(j)) = d30(index{i}(j)) + ...
                               (nParticle{i}(j) * (d{i}(j)^3));
            
            d20(index{i}(j)) = d20(index{i}(j)) + ...
                               (nParticle{i}(j) * (d{i}(j)^2));
            
            d10{i}(index{i}(j)) = d10{i}(index{i}(j)) + ...
                                  (nParticle{i}(j) * d{i}(j));
        end
        
        % Calculate Derived Variables
        density{i} = (1000 * density{i}) / cellVolume;
        d32{i} = (d30 ./ d20) * 1e6;
        d10{i} = (d10{i} ./ nParticles{i}) * 1e6;
        
        % Set Empty Cells Back to Zero
        d32{i}(isnan(d32{i})) = 0;
        d10{i}(isnan(d10{i})) = 0;
    end
    
    % Make Arrays Sparse
    nParticles{i} = sparse(nParticles{i});
    density{i} = sparse(density{i});
    d32{i} = sparse(d32{i});
    d10{i} = sparse(d10{i});
    
    % Remove Unnecessary Data
    index{i} = [];
    nParticle{i} = [];
    d{i} = [];

    % Update Waitbar
    send(dQ, []);
end
clear index nParticle d d_tmp;

delete(wB);

volumeData.inst.nParticles = nParticles; clear nParticles;
volumeData.inst.density = density; clear density;
volumeData.inst.d32 = d32; clear d32;
volumeData.inst.d10 = d10; clear d10;

disp(' ');

% Generate Time-Averaged Volume Field
disp('    Generating Time-Averaged Volume Field...');

% Calculate Instantaneous Field Variables
disp('        Calculating Time-Averaged Field Variables');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Time-Averaged Field Variables', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Perform Calculation
volumeData.mean.nParticles = zeros([nCells,1]);
volumeData.mean.density = volumeData.mean.nParticles;
volumeData.mean.d32 = volumeData.mean.nParticles;
volumeData.mean.d10 = volumeData.mean.nParticles;

for i = 1:nTimes
    volumeData.mean.nParticles = volumeData.mean.nParticles + volumeData.inst.nParticles{i};
    volumeData.mean.density = volumeData.mean.density + volumeData.inst.density{i};
    volumeData.mean.d32 = volumeData.mean.d32 + volumeData.inst.d32{i};
    volumeData.mean.d10 = volumeData.mean.d10 + volumeData.inst.d10{i};

    % Update Waitbar
    waitbar((i / nTimes), wB);
end
clear i;

volumeData.mean.nParticles = sparse(volumeData.mean.nParticles / nTimes);
volumeData.mean.density = sparse(volumeData.mean.density / nTimes);
volumeData.mean.d32 = sparse(volumeData.mean.d32 / nTimes);
volumeData.mean.d10 = sparse(volumeData.mean.d10 / nTimes);

delete(wB);

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
    selection = input('Plot Time-Averaged Volume Field? [y/n]: ', 's');

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
    selection = input('Plot Instantaneous Volume Field? [y/n]: ', 's');

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


%% Present Volume Fields

disp('Volume Field Presentation');
disp('--------------------------');

disp(' ');

if plotInst || plotMean
    gridShape = [height(unique(volumeData.positionGrid(:,1))), ...
                 height(unique(volumeData.positionGrid(:,2))), ...
                 height(unique(volumeData.positionGrid(:,3)))];
             
    spatialRes = cellSize.target / 2;
    xInit = reshape(volumeData.positionGrid(:,1), gridShape);
    yInit = reshape(volumeData.positionGrid(:,2), gridShape);
    zInit = reshape(volumeData.positionGrid(:,3), gridShape);
    POD = false;
    
    if contains(caseID, 'Windsor')
        cMap = viridis(3);
        
        if strcmp(caseID, 'Windsor_SB_wW_Upstream_SC') || strcmp(caseID, 'Windsor_SB_fullScale_multiPhase')
            cMap = cMap(1,:);
        elseif strcmp(caseID, 'Windsor_ST_wW_Upstream_SC')
            cMap = cMap(2,:);
        elseif strcmp(caseID, 'Windsor_RSST_wW_Upstream_SC')
            cMap = cMap(3,:);
        end
        
    else
        cMap = viridis(1);
    end

    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    
    switch format
        
        case 'A' % 0.75 L
            
            if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
                xLimsPlot = [0.31875; 1.43075];
                yLimsPlot = [-0.4176; 0.4176];
                zLimsPlot = [0; 0.4176];
            else
                xLimsPlot = [1.275; 5.723];
                yLimsPlot = [-1.6704; 1.6704];
                zLimsPlot = [0; 1.6704];
            end
            
        case 'B' % 2 L
            
            if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
                xLimsPlot = [0.31875; 2.73575];
                yLimsPlot = [-0.522; 0.522];
                zLimsPlot = [0; 0.522];
            else
                xLimsPlot = [1.275; 10.943];
                yLimsPlot = [-2.088; 2.088];
                zLimsPlot = [0; 2.088];
            end
            
        case 'C' % 4 L
            
            xLimsPlot = [1.275; 19.295];
            yLimsPlot = [-2.5056; 2.5056];
            zLimsPlot = [0; 2.5056];
            
    end
    
    if normalise
        xLimsPlot = round((xLimsPlot / normLength), spacePrecision);
        yLimsPlot = round((yLimsPlot / normLength), spacePrecision);
        zLimsPlot = round((zLimsPlot / normLength), spacePrecision);
    end
    
end

if plotMean
    disp('    Presenting Time-Averaged Volume Field...');
    
    fieldData = reshape(full(volumeData.mean.density), gridShape);
    
    if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
        isoValue = [];
    else
        isoValue = [5e-2; 5e-4];
    end
    
    figSubtitle = ' ';
    
    viewAngle = [0, 0;
                 0, 90;
                 90, 0;
                 30, 30];
    
    for i = 1:height(isoValue)
        
        switch format

            case 'A'
                figName = ['Near_Wake_Average_Spray_Density_', num2str(isoValue(i)), '_kg_m3'];

            case 'B'
                figName = ['Mid_Wake_Average_Spray_Density_', num2str(isoValue(i)), '_kg_m3'];

            case 'C'
                figName = ['Far_Wake_Average_Spray_Density_', num2str(isoValue(i)), '_kg_m3'];

        end

        fig = plotVolumeField(xLimsData, yLimsData, zLimsData, spatialRes, xInit, yInit, zInit, ...
                              POD, fieldData, fig, figName, geometry, isoValue(i), cMap, figTitle, ...
                              figSubtitle, viewAngle, xLimsPlot, yLimsPlot, zLimsPlot, figSave);
        
    end
    clear i;
                       
    disp(' ');
end

if plotInst
    disp('    Presenting Instantaneous Volume Field...');
    
    if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
        isoValue = [];
    else
        isoValue = [1e0; 1e-1];
    end
    
    viewAngle = [30, 30];
    
    for i = 1:height(isoValue)
        figHold = fig;
    
        for j = startFrame:endFrame

            if j ~= startFrame
                clf(fig);
                fig = figHold;
            end
            
            fieldData = reshape(full(volumeData.inst.density{j}), gridShape);
            figTime = num2str(volumeData.time(j), ['%.', num2str(timePrecision), 'f']);

            switch format

                case 'A'
                    figName = ['Near_Wake_Inst_Spray_Density_', num2str(isoValue(i)), ...
                               'kg_m3_T', erase(figTime, '.')];
                    
                case 'B'
                    figName = ['Mid_Wake_Inst_Spray_Density_', num2str(isoValue(i)), ...
                               'kg_m3_T', erase(figTime, '.')];
                    
                case 'C'
                    figName = ['Far_Wake_Inst_Spray_Density_', num2str(isoValue(i)), ...
                               'kg_m3_T', erase(figTime, '.')];

            end

            figSubtitle = [figTime, ' \it{s}'];

            fig = plotVolumeField(xLimsData, yLimsData, zLimsData, spatialRes, xInit, yInit, zInit, ...
                                  POD, fieldData, fig, figName, geometry, isoValue(i), cMap, figTitle, ...
                                  figSubtitle, viewAngle, xLimsPlot, yLimsPlot, zLimsPlot, figSave);

        end
        clear j;
        
    end
    clear i;
    
    disp(' ');
end

if ~plotMean && ~plotInst
    disp('    Skipping Volume Field Presentation');

    disp(' ');
end

disp(' ');


%% Save Volume Field Data

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
                
                if ~exist([saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/nearWake'], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/nearWake']);
                end
                
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/nearWake/', dataID, '.mat']);
                
                save([saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/nearWake/', dataID, '.mat'], ...
                      'caseID', 'dataID', 'volumeData', 'cellSize', 'sampleInterval', 'timePrecision', 'dLims', 'normalise', '-v7.3', '-noCompression');
                
                disp('        Success');
                
            case 'B'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/midWake'], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/midWake']);
                end
                
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/midWake/', dataID, '.mat']);
                
                save([saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/midWake/', dataID, '.mat'], ...
                      'caseID', 'dataID', 'volumeData', 'cellSize', 'sampleInterval', 'timePrecision', 'dLims', 'normalise', '-v7.3', '-noCompression');
                
                disp('        Success');
                
            case 'C'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/farWake'], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/farWake']);
                end
                
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/farWake/', dataID, '.mat']);
                
                save([saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/farWake/', dataID, '.mat'], ...
                      'caseID', 'dataID', 'volumeData', 'cellSize', 'sampleInterval', 'timePrecision', 'dLims', 'normalise', '-v7.3', '-noCompression');
                
                disp('        Success');
                
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;


%% Local Functions

function D = inputD(type)

    D = str2double(input(['    ', type, ' Diameter of Interest [', char(956), 'm]: '], 's'));
    
    if isnan(D) || length(D) > 1 || D < 1
        disp('        WARNING: Invalid Entry');
        
        D = -1;
    end
    
end


function frameNo = inputFrames(Nt, type)

    frameNo = str2double(input(['    Input Desired ', type, ' Frame [1-', num2str(Nt), ']: '], 's'));
    
    if isnan(frameNo) || frameNo < 1 || frameNo > Nt
        disp('        WARNING: Invalid Entry');
        
        frameNo = -1;
    end

end        
