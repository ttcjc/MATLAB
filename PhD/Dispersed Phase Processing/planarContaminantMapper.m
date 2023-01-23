%% Planar Lagrangian Contaminant Mapper v2.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

saveLocation = '/mnt/Processing/Data';
% saveLocation = '~/Data';

normalise = true; % Normalisation of Dimensions

cloudName = 'kinematicCloud'; % OpenFOAM Cloud Name

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

cellSize = 8e-3; % Spatial Resolution of Contaminant Map [m or l]

normalisationValue = 7.021451812671253e-09; % Windsor_SB_wW_Upstream_SC 1L Time-Averaged Max

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp('==============================');
disp('Planar Contaminant Mapper v2.0');
disp('==============================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit (Base Contamination and Far-Field Extraction Plane)
% v1.1 - Updated Name and Added Time-Averaging Functionality
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats


%% Initialise Case

[caseFolder, caseName, timeDirs, deltaT, timePrecision, geometry, ...
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
        [dataID, LagProps, ~, LagData, ~, sampleInterval] = initialiseLagData(saveLocation, caseFolder, caseName, ...
                                                                              cloudName, false, true, ...
                                                                              false, timeDirs, deltaT, ...
                                                                              timePrecision, nProc);
                                                                                
    case 'B'
        [dataID, LagProps, LagData, ~, ~, sampleInterval] = initialiseLagData(saveLocation, caseFolder, caseName, ...
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
        planePos = erase(planes{index}, '.');
        clear planes;
        
        disp(['Plane of Interest: ', planePos]);
end

if normalise
    dataID = [dataID, '_Norm'];
end

disp(' ');
disp(' ');


%% Select Mapping Options

disp('Mapping Options');
disp('----------------');

dLims = zeros(2,1);

valid = false;
while ~valid
    disp(' ');
    selection = input('Filter Particle Diameters? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        dLims = [1; 120];
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
        
        if (dLims(2) < 1) || (dLims(1) > 120)
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
    dataID = insertBefore(dataID, '_Norm', ['_D', num2str(dLims(1)), '_D', num2str(dLims(2))]);
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
evalc('parpool(nProc);');

disp(' ');

disp('    Initialising...');

% Identify Empty Time Instances
i = 1;
while i <= height(LagData.time)
    
    if isempty(LagData.timeExact{i})
        LagData.timeExact{i} = -1;
        
        for j = 1:height(LagProps)
            LagData.(LagProps{j}){i} = -1;
        end
        
    else
        i = i + 1;
    end
    
end
clear i;

% Shift Data Origin
if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
    
    for i = 1:height(LagData.time)
        
        if LagData.positionCartesian{i} ~= -1
            LagData.positionCartesian{i}(:,1) = LagData.positionCartesian{i}(:,1) + 1.325;
        end
        
    end
    
end

% Normalise Dimensions
if normalise
    
    if contains(caseName, ["Run_Test", "Windsor"])
        
        for i = 1:height(LagData.time)
            
            if LagData.positionCartesian{i} ~= -1
                LagData.positionCartesian{i}  = round((LagData.positionCartesian{i} / 1.044), ...
                                                      spacePrecision);
            end
            
        end
        
    end
    
end

% Specify Map Boundaries
switch format
    
    case 'A'
        % Identify Model Base
        parts = fieldnames(geometry);
        for i = 1:height(parts)
            
            if max(geometry.(parts{i}).vertices(:,1)) == xDims(2)
                break;
            end
            
            if i == height(parts)
                error('Mismatch Between ''xDims'' and Geometry Bounding Box')
            end
            
        end
    
        geoPoints = geometry.(parts{i}).vertices;
        basePoints = geoPoints((geoPoints(:,1) == xDims(2)),:);
        
        mapPerim = boundary(basePoints(:,2), basePoints(:,3), 0.95);
        mapPerim = basePoints(mapPerim,:);
        basePoly = polyshape(mapPerim(:,2), mapPerim(:,3), 'keepCollinearPoints', true);
        basePoly = polybuffer(basePoly, -0.0025, 'jointType', 'square');
        mapPerim = ones(height(basePoly.Vertices),3) * mapPerim(1,1);
        mapPerim(:,[2,3]) = basePoly.Vertices(:,[1,2]);

        if ~all(mapPerim(1,:) == mapPerim(end,:))
            mapPerim = [mapPerim; mapPerim(1,:)]; % Close Boundary
        end
        
        clear basePoints basePoly;
        
        xLimsData = xDims(2);
        yLimsData = [min(mapPerim(:,2)); max(mapPerim(:,2))];
        zLimsData = [min(mapPerim(:,3)); max(mapPerim(:,3))];
        
    case 'B'
    
        if contains(caseName, ["Run_Test", "Windsor"])
            mapPerim = [];
            
            xLimsData = LagData.positionCartesian{end}(1,1);
            yLimsData = [-0.5945; 0.5945];
            zLimsData = [0; 0.739];
            
            if normalise
                yLimsData = round((yLimsData / 1.044), spacePrecision);
                zLimsData = round((zLimsData / 1.044), spacePrecision);
            end
            
        end
        
end

disp(' ');

% Identify Particles of Interest
disp('    Identifying Particles of Interest...');

% Initialise Progress Bar
wB = waitbar(0, 'Collating Particles of Interest', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, height(LagData.time));

% Collate Particles of Interest
index = cell(height(LagData.time),1);

d = LagData.d;
positionCartesian = LagData.positionCartesian;
parfor i = 1:height(LagData.time)
    
    if positionCartesian{i} ~= -1
        index{i} = find(((d{i} * 1e6) >= dLims(1)) & ...
                        ((d{i} * 1e6) <= dLims(2)) & ...
                        (positionCartesian{i}(:,1) == xLimsData) & ...
                        (positionCartesian{i}(:,2) >= yLimsData (1)) & ...
                        (positionCartesian{i}(:,2) <= yLimsData (2)) & ...
                        (positionCartesian{i}(:,3) >= zLimsData (1)) & ...
                        (positionCartesian{i}(:,3) <= zLimsData (2))); %#ok<PFBNS>
    end
    
    send(dQ, []);
end
clear d positionCartesian;

delete(wB);

% Remove Unnecessary Data
for i = 1:height(LagData.time)
    LagData.timeExact{i} = LagData.timeExact{i}(index{i});
    
    for j = 1:height(LagProps)
        LagData.(LagProps{j}){i} = LagData.(LagProps{j}){i}(index{i},:);
    end

end

disp(' ');

% Generate Instantaneous Contaminant Maps
disp('    Generating Instantaneous Contaminant Maps...');

% Adjust Uniform Cell Size to Fit Region of Interest
cellSizeX = cellSize;
cellSizeY = (yLimsData (2) - yLimsData (1)) / round(((yLimsData (2) - yLimsData (1)) / cellSize));
cellSizeZ = (zLimsData (2) - zLimsData (1)) / round(((zLimsData (2) - zLimsData (1)) / cellSize));

[y, z] = ndgrid(yLimsData(1):cellSizeY:yLimsData(2), zLimsData(1):cellSizeZ:zLimsData(2));

mapData.positionGrid = zeros(height(y(:)),3);
mapData.positionGrid(:,1) = xLimsData;
mapData.positionGrid(:,(2:3)) = [y(:), z(:)];

mapData.inst.time = LagData.time;

% Initialise Progress Bar
wB = waitbar(0, 'Assigning Particles to Map Nodes', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, height(mapData.inst.time));

% Assign Particles to Map Nodes
index = cell(height(mapData.inst.time),1); % Array Position of Closest Mesh Node

positionGrid = mapData.positionGrid;
positionCartesian = LagData.positionCartesian;
parfor i = 1:height(mapData.inst.time)
    
    if positionCartesian{i} ~= -1
        index{i} = dsearchn(positionGrid, positionCartesian{i});
    end
    
    send(dQ, []);
end
clear positionGrid positionCartesian;

delete(wB);

% Offset Particles From Base for Better Visibility
switch format
    
    case 'A'
        mapData.positionGrid(:,1) = mapData.positionGrid(:,1) + 1e-3;
        
end

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Mapping Variables', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, height(mapData.inst.time));

% Calculate Instantaneous Mapping Variables
nParticles = cell(height(LagData.time),1); % Number of Particles in Cell
d10 = nParticles; % Arithmetic Mean Diameter in Cell
d20 = nParticles; % Surface Mean Diameter in Cell
d30 = nParticles; % Volume Mean Diameter in Cell
d32 = nParticles; % Sauter Mean Diameter in Cell
mass = nParticles; % Total Mass in Cell
massNorm = nParticles; % Normalised Mass in Cell

positionGrid = mapData.positionGrid;
positionCartesian = LagData.positionCartesian;
nParticle = LagData.nParticle;
d = LagData.d;
parfor i = 1:height(mapData.inst.time)
    nParticles{i} = zeros(height(positionGrid),1);
    d10{i} = nParticles{i};
    d20{i} = nParticles{i};
    d30{i} = nParticles{i};
    d32{i} = nParticles{i};
    mass{i} = nParticles{i};
    massNorm{i} = nParticles{i};
    
    for j = 1:height(positionCartesian{i})
        nParticles{i}(index{i}(j)) = nParticles{i}(index{i}(j)) + ...
                                     nParticle{i}(j);
        d10{i}(index{i}(j)) = d10{i}(index{i}(j)) + ...
                              (nParticle{i}(j) * d{i}(j));
        d20{i}(index{i}(j)) = d20{i}(index{i}(j)) + ...
                              (nParticle{i}(j) * (d{i}(j)^2));
        d30{i}(index{i}(j)) = d30{i}(index{i}(j)) + ...
                              (nParticle{i}(j) * (d{i}(j)^3));
        mass{i}(index{i}(j)) = mass{i}(index{i}(j)) + ...
                               (nParticle{i}(j) * ((1 / 12) * tau * (d{i}(j)^3)));
    end
    
    mass{i} = 1000 * mass{i};
    massNorm{i} = mass{i} / normalisationValue;
    d32{i} = (d30{i} ./ d20{i}) * 1e6;
    d30{i} = ((d30{i} ./ nParticles{i}).^(1/3)) * 1e6;
    d20{i} = ((d20{i} ./ nParticles{i}).^(1/2)) * 1e6;
    d10{i} = (d10{i} ./ nParticles{i}) * 1e6;
    
    % Set Empty Cells Back to Zero
    d10{i}(isnan(d10{i})) = 0;
    d20{i}(isnan(d20{i})) = 0;
    d30{i}(isnan(d30{i})) = 0;
    d32{i}(isnan(d32{i})) = 0;
    
    send(dQ, []);
end
clear positionGrid positionCartesian nParticle d;

delete(wB);

mapData.inst.nParticles = nParticles;
mapData.inst.d10 = d10;
mapData.inst.d20 = d20;
mapData.inst.d30 = d30;
mapData.inst.d32 = d32;
mapData.inst.mass = mass;
mapData.inst.massNorm = massNorm;
clear nParticles d10 d20 d30 d32 mass massNorm;

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Centre of Mass', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Calculate Instantaneous Centre of Mass
mapData.inst.CoM = cell(height(mapData.inst.time),1);

for i = 1:height(mapData.inst.time)
    mapData.inst.CoM{i} = zeros(1,3);
    
    mapData.inst.CoM{i}(1) = mapData.positionGrid(1,1);
    mapData.inst.CoM{i}(2) = sum(mapData.inst.mass{i} .* mapData.positionGrid(:,2)) / sum(mapData.inst.mass{i});
    mapData.inst.CoM{i}(3) = sum(mapData.inst.mass{i} .* mapData.positionGrid(:,3)) / sum(mapData.inst.mass{i});
    
    waitbar((i / height(mapData.inst.time)), wB);
end

delete(wB);

disp(' ');

% Generate Time-Averaged Contaminant Map
disp('    Generating Time-Averaged Contaminant Maps...');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Time-Averaged Mapping Variables', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, height(mapData.inst.time));

% Calculate Time-Averaged Mapping Variables
nParticlesMean = zeros(height(mapData.positionGrid),1);
d10Mean = nParticlesMean;
d20Mean = nParticlesMean;
d30Mean = nParticlesMean;
d32Mean = nParticlesMean;
massMean = nParticlesMean;

nParticles = mapData.inst.nParticles;
d10 = mapData.inst.d10;
d20 = mapData.inst.d20;
d30 = mapData.inst.d30;
d32 = mapData.inst.d32;
mass = mapData.inst.mass;
parfor i = 1:height(mapData.inst.time)
    nParticlesMean = nParticlesMean + nParticles{i};
    d10Mean = d10Mean + d10{i};
    d20Mean = d20Mean + d20{i};
    d30Mean = d30Mean + d30{i};
    d32Mean = d32Mean + d32{i};
    massMean = massMean + mass{i};
    
    send(dQ, []);
end
clear nParticles d10 d20 d30 d32 mass;

delete(wB);

mapData.mean.nParticles = nParticlesMean / height(mapData.inst.time);
mapData.mean.d10 = d10Mean / height(mapData.inst.time);
mapData.mean.d20 = d20Mean / height(mapData.inst.time);
mapData.mean.d30 = d30Mean / height(mapData.inst.time);
mapData.mean.d32 = d32Mean / height(mapData.inst.time);
mapData.mean.mass = massMean / height(mapData.inst.time);
mapData.mean.massNorm = mapData.mean.mass / normalisationValue;
clear nParticlesMean d10Mean d20Mean d30Mean d32Mean massMean;

% Calculate Time-Averaged Centre of Mass
mapData.mean.CoM = zeros(1,3);
    
mapData.mean.CoM(1) = mapData.positionGrid(1,1);
mapData.mean.CoM(2) = sum(mapData.mean.mass .* mapData.positionGrid(:,2)) / sum(mapData.mean.mass);
mapData.mean.CoM(3) = sum(mapData.mean.mass .* mapData.positionGrid(:,3)) / sum(mapData.mean.mass);

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
        nFrames = inputFrames(height(mapData.inst.time));
        
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

            if contains(caseName, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 1.52725];
                yLimsPlot = [-0.2445; 0.2445];
                zLimsPlot = [0; 0.389];
            end

        case 'B'
            orientation = 'YZ';

            if contains(caseName, ["Run_Test", "Windsor"])
%                 % Plot Numerical Data Range
%                 xLimsPlot = [0.31875; 4.65925];
%                 yLimsPlot = [-0.5945; 0.5945];
%                 zLimsPlot = [0; 0.739];
                
                % Plot Experimental Data Range
                xLimsPlot = [0.31875; 4.65925];
                yLimsPlot = [-0.399; 0.218];
                zLimsPlot = [0.0105; 0.4985];
            end

    end

    if normalise
        xLimsPlot = round((xLimsPlot / 1.044), spacePrecision);
        yLimsPlot = round((yLimsPlot / 1.044), spacePrecision);
        zLimsPlot = round((zLimsPlot / 1.044), spacePrecision);
    end

    positionData = mapData.positionGrid;
    cMap = flipud(viridis(32));
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
end

if plotMean
    
    for i = 1:height(plotVars)
        disp(['    Presenting Time-Averaged ''', plotVars{i}, ''' Data...']);
        
        scalarData = mapData.mean.(plotVars{i});
        
        switch format
            
            case 'A'
                figName = ['Time_Averaged_Base_', plotVars{i}, '_Map', '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];
                
            case 'B'
                figName = ['Time_Averaged_', planePos, '_', plotVars{i}, '_Map', '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];
                
        end
        
        contourlines = (0.2:0.2:0.8);
        
        if any(strcmp(plotVars{i}, {'mass', 'massNorm'}))
            CoM = mapData.mean.CoM;
        else
            CoM = [];
        end
            
        figSubtitle = ' ';
        
        if any(strcmp(plotVars{i}, {'d10', 'd20', 'd30', 'd32'}))
%             cLims = dLims;
            cLims = [0; 40]; % Max Planar Contamination
        elseif strcmp(plotVars{i}, 'massNorm')
            cLims = [0; 1]; % Max Base Contamination
%             cLims = [0; 20]; % Max Planar Contamination
        else
            cLims = [0; max(scalarData)];
        end
        
        fig = planarScalarPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                mapPerim, fig, figName, cMap, geometry, contourlines, ...
                                xDims, yDims, zDims, CoM, figTitle, figSubtitle, cLims, ...
                                xLimsPlot, yLimsPlot, zLimsPlot, normalise);
    end
    
    disp(' ');
end

if plotInst
    
    for i = 1:height(plotVars)
        disp(['    Presenting Instantaneous ''', plotVars{i}, ''' Data...']);
        
        contourlines = [];
        
        if any(strcmp(plotVars{i}, {'d10', 'd20', 'd30', 'd32'}))
            cLims = dLims;
        elseif strcmp(plotVars{i}, 'massNorm')
            cLims = [0; 3.6];
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
            figTime = num2str(mapData.inst.time(j), ['%.', num2str(timePrecision), 'f']);
            
            switch format
                
                case 'A'
                    figName = ['Instantaneous_Base_', plotVars{i}, '_Map_T', erase(figTime, '.'), '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];
                
                case 'B'
                    figName = ['Instantaneous_', planePos, '_', plotVars{i}, '_Map_T', erase(figTime, '.'), '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];
            end
            
        if contains(plotVars{i}, ["mass", "massNorm"])
            CoM = mapData.inst.CoM{j};
        else
            CoM = [];
        end
            
        figSubtitle = [figTime, ' \it{s}'];
        
        fig = planarScalarPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                mapPerim, fig, figName, cMap, geometry, contourlines, ...
                                xDims, yDims, zDims, CoM, figTitle, figSubtitle, cLims, ...
                                xLimsPlot, yLimsPlot, zLimsPlot, normalise);
        end
        
    end
    
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
                
                if ~exist([saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseName, '/base'], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseName, '/base']);
                end
                
            case 'B'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseName, '/', planePos], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseName, '/', planePos]);
                end
                
        end
        
        switch format
            
            case 'A'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseName, '/base/', dataID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseName, '/base/', dataID, '.mat'], ...
                     'dataID', 'mapData', 'sampleInterval', 'dLims', 'normalise', 'timePrecision', '-v7.3', '-noCompression');
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseName, '/', planePos, '/', dataID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseName, '/', planePos, '/', dataID, '.mat'], ...
                     'dataID', 'mapData', 'sampleInterval', 'dLims', 'normalise', 'timePrecision', '-v7.3', '-noCompression');
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