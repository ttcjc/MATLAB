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


%% Planar Lagrangian Contaminant Mapper v2.3

cellSize.target = 8e-3; % Spatial Resolution of Contaminant Map [m or l]

cloudName = 'kinematicCloud'; % OpenFOAM Cloud Name

figSave = false; % Save .fig File(s)

normalise = true; % Normalisation of Dimensions

% normalisationValue = 2.214802385935855e-08; % Value Used to Normalise Masses (SB, 0.5 L, 200 Hz)
normalisationValue = 1.1431747e-08; % Value Used to Normalise Masses (SB, 1.0 L, 200 Hz)

disp('==============================');
disp('Planar Contaminant Mapper v2.3');
disp('==============================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit (Base Contamination and Far-Field Extraction Plane)
% v1.1 - Updated Name and Added Time-Averaging Functionality
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats
% v2.1 - Improved Efficiency of Map Generation
% v2.2 - Added Support for Multiple Mean Particle Diameter Definitions (Disabled by Default)
% v2.3 - Shifted to Thread-Based Parallelization to Reduce Memory Requirements


%% Initialise Case

[caseFolder, caseID, timeDirs, deltaT, timePrecision, geometry, ...
 xDims, yDims, zDims, spacePrecision, normalise] = initialiseCaseData(normalise);

disp(' ');
disp(' ');


%% Select Mapping Location

disp('Mapping Location');
disp('-----------------');

disp(' ');

disp('Possible Mapping Locations:');
disp('    A: Base Contamination');
disp('    B: Far-Field Spray Transport');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Mapping Location [A/B]: ', 's');

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


%% Initialise Lagrangian Data

switch format
    
    case 'A'
        [dataID, LagProps, ~, LagData, ~, sampleInterval] = initialiseLagData(saveLocation, caseFolder, caseID, ...
                                                                              cloudName, false, true, ...
                                                                              false, timeDirs, deltaT, ...
                                                                              timePrecision, nProc);
                                                                                
    case 'B'
        [dataID, LagProps, LagData, ~, ~, sampleInterval] = initialiseLagData(saveLocation, caseFolder, caseID, ...
                                                                              cloudName, true, false, ...
                                                                              false, timeDirs, deltaT, ...
                                                                              timePrecision, nProc);

        % Select Plane of Interest
        planes = fieldnames(LagData);
        
        valid = false;
        while ~valid
            disp(' ');
            [index, valid] = listdlg('listSize', [300, 300], ...
                                     'selectionMode', 'single', ...
                                     'name', 'Select Plane of Interest', ...
                                     'listString', planes);
            
            if ~valid
                disp('WARNING: No Plane of Interest Selected');
            end
        
        end
        clear valid;

        LagData = LagData.(planes{index});
        planeID = erase(planes{index}, '.');
        
        clear planes;
        
        disp(['Plane of Interest: ', planeID]);
end

disp(' ');
disp(' ');


%% Select Mapping Options

disp('Mapping Options');
disp('----------------');

dLims = zeros([2,1]);

valid = false;
while ~valid
    disp(' ');
    selection = input('Filter Particle Diameters? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        
        if contains(caseID, ["Run_Test", "Windsor"])
            dLims = [1; 147];
        end
        
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
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
        
        if contains(caseID, "Run_Test", "Windsor") && ((dLims(2) < 1) || (dLims(1) > 147))
            disp('        WARNING: No Lagrangian Data in Diameter Range');
            continue;
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

if normalise
    dataID = [dataID, '_D', num2str(dLims(1)), '_D', num2str(dLims(2)), '_Norm'];
else
    dataID = [dataID, '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];
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

if normalise
    
    % Normalise Coordinate System
    disp('        Normalising Coordinate System');
    
    if contains(caseID, ["Run_Test", "Windsor"])
        
        % Initialise Progress Bar
        wB = waitbar(0, 'Normalising Coordinate System', 'name', 'Progress');
        wB.Children.Title.Interpreter = 'none';
        
        % Perform Normalisation
        for i = 1:nTimes
            
            if ~isempty(LagData.positionCartesian{i})
                LagData.positionCartesian{i}  = round((LagData.positionCartesian{i} / 1.044), ...
                                                      spacePrecision);
            end
            
            % Update Waitbar
            waitbar((i / nTimes), wB);
        end
        clear i;
        
        delete(wB);
    else
        disp('            WARNING: Dimension Normalisation for This Case Type Not Supported');
    end
    
end

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
        
    case 'B'
    
        if contains(caseID, ["Run_Test", "Windsor"])
            mapPerim = [];
            
            xLimsData = double(LagData.positionCartesian{end}(1,1));
            yLimsData = [-0.522; 0.522];
            zLimsData = [0; 0.6264];
            
            if normalise
                yLimsData = round((yLimsData / 1.044), spacePrecision);
                zLimsData = round((zLimsData / 1.044), spacePrecision);
            end
            
        end
        
end

disp(' ');

% Collate Particles of Interest
disp('    Collating Particles of Interest...');

index = cell(nTimes,1);

switch format

    case 'A'
        
        for i = 1:nTimes
            
            if ~isempty(LagData.positionCartesian{i})
                index{i} = find(((LagData.d{i} * 1e6) >= dLims(1)) & ...
                                ((LagData.d{i} * 1e6) <= dLims(2)) & ...
                                (round(LagData.positionCartesian{i}(:,1), spacePrecision) == single(round(xLimsData, spacePrecision))) & ...
                                (LagData.positionCartesian{i}(:,2) >= yLimsData(1)) & ...
                                (LagData.positionCartesian{i}(:,2) <= yLimsData(2)) & ...
                                (LagData.positionCartesian{i}(:,3) >= zLimsData(1)) & ...
                                (LagData.positionCartesian{i}(:,3) <= zLimsData(2)));
            end
        
        end
        clear i;

    case 'B'
        
        for i = 1:nTimes
            
            if ~isempty(LagData.positionCartesian{i})
                index{i} = find(((LagData.d{i} * 1e6) >= dLims(1)) & ...
                                ((LagData.d{i} * 1e6) <= dLims(2)) & ...
                                (LagData.positionCartesian{i}(:,2) >= yLimsData(1)) & ...
                                (LagData.positionCartesian{i}(:,2) <= yLimsData(2)) & ...
                                (LagData.positionCartesian{i}(:,3) >= zLimsData(1)) & ...
                                (LagData.positionCartesian{i}(:,3) <= zLimsData(2)));
            end
        
        end
        clear i;

end

% Remove Unnecessary Data
disp('        Removing Unnecessary Data');

LagFields = fieldnames(LagData);
reqFields = {'time'; 'timeExact'; 'd'; 'nParticle'; 'positionCartesian'};

LagData = rmfield(LagData, LagFields(~ismember(LagFields, reqFields)));
LagProps = LagProps(ismember(LagProps, reqFields));

for i = 1:nTimes
    LagData.timeExact{i} = LagData.timeExact{i}(index{i});
    
    for j = 1:height(LagProps)
        LagData.(LagProps{j}){i} = LagData.(LagProps{j}){i}(index{i},:);
    end
    clear j;

end
clear i;

disp(' ');

% Generate Presentation Grid
disp('    Generating 2D Presentation Grid...');

% Adjust Uniform Cell Size to Fit Region of Interest
if normalise
    cellSize.target = round((cellSize.target / 1.044), spacePrecision);
end

cellSize.x = cellSize.target;
cellSize.y = (yLimsData (2) - yLimsData (1)) / round(((yLimsData (2) - yLimsData (1)) / cellSize.target));
cellSize.z = (zLimsData (2) - zLimsData (1)) / round(((zLimsData (2) - zLimsData (1)) / cellSize.target));

cellSize.area = cellSize.y * cellSize.z;

[y, z] = ndgrid(yLimsData(1):cellSize.y:yLimsData(2), zLimsData(1):cellSize.z:zLimsData(2));

mapData.positionGrid = zeros([height(y(:)),3]);
mapData.positionGrid(:,1) = xLimsData;
mapData.positionGrid(:,(2:3)) = [y(:), z(:)];

% Assign Particles to Grid Cells
disp('        Assigning Particles to Nearest Grid Cell');

% Initialise Progress Bar
wB = waitbar(0, 'Assigning Particles to Nearest Grid Cell', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, nTimes);

% Perform Assignment
totalParticles = cellfun(@height, LagData.positionCartesian);
index = cell(nTimes,1); % Array Position of Closest Mesh Node

positionGrid = mapData.positionGrid;
positionCartesian = LagData.positionCartesian;
parfor i = 1:nTimes
    
    if totalParticles(i) > 0
        index{i} = dsearchn(positionGrid, positionCartesian{i});
    end

    % Remove Unnecessary Data
    positionCartesian{i} = [];

    % Update Waitbar
    send(dQ, []);
end
clear i positionGrid positionCartesian;

delete(wB);

disp(' ');

% Generate Instantaneous Contaminant Maps
disp('    Generating Instantaneous Contaminant Maps...');

mapData.time = LagData.time;

% Calculate Instantaneous Field Variables
disp('        Calculating Instantaneous Field Variables');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Field Variables', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, nTimes);

% Perform Calculation
nParticles = cell(nTimes,1); nParticles(:) = {zeros([height(mapData.positionGrid),1], 'single')}; % Number of Particles in Cell
mass = nParticles; % Total Mass in Cell
massNorm = nParticles; % Normalised Mass in Cell
% d32 = nParticles; % Sauter Mean Diameter in Cell
% d30 = nParticles; % Volume Mean Diameter in Cell
d20 = nParticles; % Surface Mean Diameter in Cell
d10 = nParticles; % Arithmetic Mean Diameter in Cell

nParticle = LagData.nParticle;
d = LagData.d;
parfor i = 1:nTimes
    
    if totalParticles(i) > 0
        
        for j = 1:totalParticles(i)
            nParticles{i}(index{i}(j)) = nParticles{i}(index{i}(j)) + ...
                                         nParticle{i}(j);
            
            mass{i}(index{i}(j)) = mass{i}(index{i}(j)) + ...
                                   (nParticle{i}(j) * ((1 / 12) * tau * (d{i}(j)^3)));
            
%             d30{i}(index{i}(j)) = d30{i}(index{i}(j)) + ...
%                                   (nParticle{i}(j) * (d{i}(j)^3));
            
            d20{i}(index{i}(j)) = d20{i}(index{i}(j)) + ...
                                  (nParticle{i}(j) * (d{i}(j)^2));
            
            d10{i}(index{i}(j)) = d10{i}(index{i}(j)) + ...
                                  (nParticle{i}(j) * d{i}(j));
        end
        
    end

    % Remove Unnecessary Data
    index{i} = [];
    nParticle{i} = [];
    d{i} = [];

    % Calculate Derived Variables
    mass{i} = 1000 * mass{i};
    massNorm{i} = mass{i} / normalisationValue;
%     d32{i} = (d30{i} / d20{i}) * 1e6;
%     d30{i} = ((d30{i} ./ nParticles{i}).^(1/3)) * 1e6;
    d20{i} = ((d20{i} ./ nParticles{i}).^(1/2)) * 1e6;
    d10{i} = (d10{i} ./ nParticles{i}) * 1e6;
    
    % Set Empty Cells Back to Zero
%     d32{i}(isnan(d32{i})) = 0;
%     d30{i}(isnan(d30{i})) = 0;
    d20{i}(isnan(d20{i})) = 0;
    d10{i}(isnan(d10{i})) = 0;
    
    % Update Waitbar
    send(dQ, []);
end
clear i j positionCartesian nParticle d;

delete(wB);

clear LagData

mapData.inst.nParticles = nParticles; clear nParticles;
mapData.inst.mass = mass; clear mass;
mapData.inst.massNorm = massNorm; clear massNorm;
% mapData.inst.d32 = d32; clear d32;
% mapData.inst.d30 = d30; clear d30;
mapData.inst.d20 = d20; clear d20;
mapData.inst.d10 = d10; clear d10;

% Calculate Instantaneous Centre of Mass
disp('        Calculating Instantaneous Centre of Mass');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Centre of Mass', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Perform Calculation
mapData.inst.CoM = cell(nTimes,1); mapData.inst.CoM(:) = {zeros([1,3], 'single')};

for i = 1:nTimes
    mapData.inst.CoM{i}(1) = mapData.positionGrid(1,1);
    mapData.inst.CoM{i}(2) = sum(mapData.inst.mass{i} .* mapData.positionGrid(:,2)) / sum(mapData.inst.mass{i});
    mapData.inst.CoM{i}(3) = sum(mapData.inst.mass{i} .* mapData.positionGrid(:,3)) / sum(mapData.inst.mass{i});
    
    % Update Waitbar
    waitbar((i / nTimes), wB);
end
clear i;

delete(wB);

disp(' ');

% Generate Time-Averaged Contaminant Map
disp('    Generating Time-Averaged Contaminant Maps...');

% Calculate Time-Averaged Field Variables
disp('        Calculating Instantaneous Field Variables');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Time-Averaged Field Variables', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Perform Calculation
mapData.mean.nParticles = zeros([height(mapData.positionGrid),1], 'single');
mapData.mean.mass = mapData.mean.nParticles;
mapData.mean.massNorm = mapData.mean.nParticles;
% mapData.mean.d32 = mapData.mean.nParticles;
% mapData.mean.d30 = mapData.mean.nParticles;
mapData.mean.d20 = mapData.mean.nParticles;
mapData.mean.d10 = mapData.mean.nParticles;

for i = 1:nTimes
    mapData.mean.nParticles = mapData.mean.nParticles + mapData.inst.nParticles{i};
    mapData.mean.mass = mapData.mean.mass + mapData.inst.mass{i};
    mapData.mean.massNorm = mapData.mean.massNorm + mapData.inst.massNorm{i};
%     mapData.mean.d32 = mapData.mean.d32 + mapData.inst.d32{i};
%     mapData.mean.d30 = mapData.mean.d30 + mapData.inst.d30{i};
    mapData.mean.d20 = mapData.mean.d20 + mapData.inst.d20{i};
    mapData.mean.d10 = mapData.mean.d10 + mapData.inst.d10{i};
    
    % Update Waitbar
    waitbar((i / nTimes), wB);
end
clear i;

delete(wB);

mapData.mean.nParticles = mapData.mean.nParticles / nTimes;
mapData.mean.mass = mapData.mean.mass / nTimes;
mapData.mean.massNorm = mapData.mean.massNorm / nTimes;
% mapData.mean.d32 = mapData.mean.d32 / nTimes;
% mapData.mean.d30 = mapData.mean.d30 / nTimes;
mapData.mean.d20 = mapData.mean.d20 / nTimes;
mapData.mean.d10 = mapData.mean.d10 / nTimes;

% Perform Calculation
disp('        Time-Averaged Centre of Mass');

mapData.mean.CoM = zeros([1,3], 'single');
    
mapData.mean.CoM(1) = mapData.positionGrid(1,1);
mapData.mean.CoM(2) = sum(mapData.mean.mass .* mapData.positionGrid(:,2)) / sum(mapData.mean.mass);
mapData.mean.CoM(3) = sum(mapData.mean.mass .* mapData.positionGrid(:,3)) / sum(mapData.mean.mass);

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

if plotInst || plotMean
    
    % Select Variable(s) of Interest
    plotVars = fieldnames(mapData.mean);
    plotVars = plotVars(1:(end - 1));

    valid = false;
    while ~valid
        disp(' ');
        [index, valid] = listdlg('listSize', [300, 300], ...
                                 'selectionMode', 'multiple', ...
                                 'name', 'Select Variable(s) to Plot', ...
                                 'listString', plotVars);

        if ~valid
            disp('WARNING: No Mapping Variables Selected');
        end
        
    end
    clear valid;

    plotVars = plotVars(index);
    
    disp(['Plotting ', num2str(height(plotVars)), ' Variable(s) of Interest']);
end

disp(' ');
disp(' ');


%% Present Contaminant Maps

disp('Map Presentation');
disp('-----------------');

disp(' ');

if plotInst || plotMean
    
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

    end
    
    if normalise
        
        if contains(caseID, ["Run_Test", "Windsor"])
            xLimsPlot = round((xLimsPlot / 1.044), spacePrecision);
            yLimsPlot = round((yLimsPlot / 1.044), spacePrecision);
            zLimsPlot = round((zLimsPlot / 1.044), spacePrecision);
        end
        
    end
    
    positionData = mapData.positionGrid;
    nPlanes = 1;
    planeNo = 1;
    cMap = flipud(viridis(32));
    refPoint = [];
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes

end

if plotMean
    
    for i = 1:height(plotVars)
        disp(['    Presenting Time-Averaged ''', plotVars{i}, ''' Data...']);
        
        scalarData = mapData.mean.(plotVars{i});
        
        switch format
            
            case 'A'
                figName = ['Time_Averaged_Base_', plotVars{i}, '_Map'];
                
            case 'B'
                figName = ['Time_Averaged_', planeID, '_', plotVars{i}, '_Map'];
                
        end
        
        if strcmp(plotVars{i}, 'massNorm')
            contourlines = [0.02; 0.02];
        else
            contourlines = [];
        end
        
        figSubtitle = ' ';
        
        if any(strcmp(plotVars{i}, {'d10', 'd20', 'd30', 'd32', 'd43'}))
            cLims = [0; 150];
        elseif strcmp(plotVars{i}, 'massNorm')
            cLims = [0; 1.2];
        else
            cLims = [0; max(scalarData)];
        end
        
        [fig, planeNo] = plotPlanarScalarField(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                               mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, contourlines, ...
                                               xDims, yDims, zDims, refPoint, figTitle, figSubtitle, cLims, ...
                                               xLimsPlot, yLimsPlot, zLimsPlot, normalise, figSave);
    end
    clear i;
    
    disp(' ');
end

if plotInst
    
    for i = 1:height(plotVars)
        disp(['    Presenting Instantaneous ''', plotVars{i}, ''' Data...']);
        
        if strcmp(plotVars{i}, 'massNorm')
            contourlines = [0.02; 0.02];
        else
            contourlines = [];
        end
        
        if any(strcmp(plotVars{i}, {'d10', 'd20', 'd30', 'd32'}))
            cLims = [0; 150];
        elseif strcmp(plotVars{i}, 'massNorm')
            cLims = [0; 4.2];
        else
            cLims = [0; max(cellfun(@max, mapData.inst.(plotVars{i})))];
        end
        
        figHold = fig;
        
        for j = 1:nFrames
            
            if j ~= 1
                clf(fig);
                fig = figHold;
            end
            
            scalarData = mapData.inst.(plotVars{i}){j};
            figTime = num2str(mapData.time(j), ['%.', num2str(timePrecision), 'f']);
            
            switch format
                
                case 'A'
                    figName = ['Instantaneous_Base_', plotVars{i}, '_Map_T', erase(figTime, '.')];
                
                case 'B'
                    figName = ['Instantaneous_', planeID, '_', plotVars{i}, '_Map_T', erase(figTime, '.')];
            end
            
            figSubtitle = [figTime, ' \it{s}'];
        
            [fig, planeNo] = plotPlanarScalarField(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                                   mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, contourlines, ...
                                                   xDims, yDims, zDims, refPoint, figTitle, figSubtitle, cLims, ...
                                                   xLimsPlot, yLimsPlot, zLimsPlot, normalise, figSave);
        end
        clear j;
        
    end
    clear i;
    
    disp(' ');
end

if ~plotMean && ~plotInst
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
        
        switch format
            
            case 'A'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseID, '/base'], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseID, '/base']);
                end
                
            case 'B'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseID, '/', planeID], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseID, '/', planeID]);
                end
                
        end
        
        switch format
            
            case 'A'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseID, '/base/', dataID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseID, '/base/', dataID, '.mat'], ...
                      'caseID', 'dataID', 'mapData', 'sampleInterval', 'timePrecision', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseID, '/', planeID, '/', dataID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseID, '/', planeID, '/', dataID, '.mat'], ...
                      'caseID', 'planeID', 'dataID', 'mapData', 'sampleInterval', 'timePrecision', 'dLims', 'normalise', '-v7.3', '-noCompression');
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


function nFrames = inputFrames(Nt)

    nFrames = str2double(input(['    Input Desired Frame Count [1-', num2str(Nt), ']: '], 's'));
    
    if isnan(nFrames) || nFrames <= 0 || nFrames > Nt
        disp('        WARNING: Invalid Entry');
        nFrames = -1;
    end

end
