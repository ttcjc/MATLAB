%% Lagrangian Contamination Mapper v2.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

normalise = true; % Normalisation of Dimensions

cloudName = 'kinematicCloud'; % OpenFOAM Cloud Name

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

cellSize = 8e-3; % Spatial Resolution of Contaminant Map [m or l]

massNormalisation = 3.810600208515075e-10; % Square-Back Base Time-Average

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp ('====================================');
disp ('Lagrangian Contamination Mapper v2.0');
disp ('====================================');

disp (' ');
disp (' ');


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
disp('    A: Surface Contamination (Base)');
disp('    B: Far-Field Spray Transport');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Mapping Location [A/B/C]: ', 's');

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
        [LagProps, ~, LagData, ~, sampleInterval] = initialiseLagData(caseFolder, caseName, cloudName, ...
                                                                      false, true, false, ...
                                                                      timeDirs, deltaT, timePrecision, nProc);
                                                                                
    case 'B'
        [LagProps, LagData, ~, ~, sampleInterval] = initialiseLagData(caseFolder, caseName, cloudName, ...
                                                                      true, false, false, ...
                                                                      timeDirs, deltaT, timePrecision, nProc);

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
                disp('WARNING: No Case Type Selected');
            end
        
        end
        clear valid;

        LagData = LagData.(planes{index});
        planePos = erase(planes{index}, '.');
        clear planes;
        
        disp(['Plane of Interest: ', planePos]);
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
            continue
        end

        dLims(2) = inputD('Max');

        if dLims(2) == -1
            continue
        end
        
        dLims = sort(dLims);
        dLims(1) = floor(dLims(1));
        dLims(2) = ceil(dLims(2));
        
        if (dLims(2) < 1) || (dLims(1) > 120)
            disp('        WARNING: No Lagrangian Data in Diameter Range');
            continue
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

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
                break
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
            mapPerim = vertcat(mapPerim, mapPerim(1,:)); % Close Boundary
        end
        
        clear basePoints basePoly;
        
        xLimsData = xDims(2);
        yLimsData = [min(mapPerim(:,2)); max(mapPerim(:,2))];
        zLimsData = [min(mapPerim(:,3)); max(mapPerim(:,3))];
        
    case 'B'
        if contains(caseName, 'Run_Test') || contains(caseName, 'Windsor')
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
                        (positionCartesian{i}(:,1) == xLimsData ) & ...
                        (positionCartesian{i}(:,2) >= yLimsData (1)) & ...
                        (positionCartesian{i}(:,2) <= yLimsData (2)) & ...
                        (positionCartesian{i}(:,3) >= zLimsData (1)) & ...
                        (positionCartesian{i}(:,3) <= zLimsData (2))); %#ok<PFBNS>
    end
    
    send(dQ, []);
end
clear d positionCartesian;

delete(wB);

contaminantData.time = LagData.time;
LagData.time = []; % Free Memory

contaminantData.timeExact = cell(height(contaminantData.time),1);

for i = 1:height(LagProps)
    contaminantData.(LagProps{i}) = contaminantData.timeExact;
end

for i = 1:height(contaminantData.time)
    contaminantData.timeExact{i} = LagData.timeExact{i}(index{i});
    LagData.timeExact{i} = []; % Free Memory
    
    for j = 1:height(LagProps)
        contaminantData.(LagProps{j}){i} = LagData.(LagProps{j}){i}(index{i},:);
        LagData.(LagProps{j}){i} = []; % Free Memory
    end
    
end

clear LagData;

disp(' ');

% Generate Instantaneous Contaminant Maps
disp('    Generating Instantaneous Contaminant Maps...');

% Adjust Uniform Cell Size to Fit Region of Interest
cellSizeX = cellSize;
cellSizeY = (yLimsData (2) - yLimsData (1)) / round(((yLimsData (2) - yLimsData (1)) / cellSize));
cellSizeZ = (zLimsData (2) - zLimsData (1)) / round(((zLimsData (2) - zLimsData (1)) / cellSize));

[y, z] = meshgrid(yLimsData (1):cellSizeY:yLimsData (2), zLimsData (1):cellSizeZ:zLimsData (2));

mapData.positionGrid = zeros(height(y(:)),3);
mapData.positionGrid(:,1) = xLimsData ;
mapData.positionGrid(:,(2:3)) = [y(:), z(:)];

mapData.inst.time = contaminantData.time;

% Initialise Progress Bar
wB = waitbar(0, 'Assigning Particles to Map Nodes', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, height(mapData.inst.time));

% Assign Particles to Map Nodes
index = cell(height(mapData.inst.time),1);

positionGrid = mapData.positionGrid;
positionCartesian = contaminantData.positionCartesian;
parfor i = 1:height(mapData.inst.time)
    
    if positionCartesian{i} ~= -1
        index{i} = dsearchn(positionGrid, positionCartesian{i});
    end
    
    send(dQ, []);
end
clear positionGrid positionCartesian;

delete(wB);

mapData.positionGrid(:,1) = mapData.positionGrid(:,1) + 1e-3;

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Mapping Variables', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, height(mapData.inst.time));

% Calculate Instantaneous Mapping Variables
nParticles = cell(height(contaminantData.time),1); % Number of Particles in Cell
d10 = nParticles; % Arithmetic Mean Diameter in Cell
d20 = nParticles; % Surface Mean Diameter in Cell
d30 = nParticles; % Volume Mean Diameter in Cell
d32 = nParticles; % Sauter Mean Diameter in Cell
mass = nParticles; % Total Mass in Cell
massNorm = nParticles; % Normalised Mass in Cell

positionGrid = mapData.positionGrid;
positionCartesian = contaminantData.positionCartesian;
nParticle = contaminantData.nParticle;
d = contaminantData.d;
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
    massNorm{i} = mass{i} / massNormalisation;
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

clear contaminantData;

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
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, height(mapData.inst.time));

% Calculate Instantaneous Centre of Mass
CoM = cell(height(mapData.inst.time),1);

mass = mapData.inst.mass;
positionGrid = mapData.positionGrid;
parfor i = 1:height(mapData.inst.time)
    CoM{i} = zeros(1,3);
    CoM{i}(1) = positionGrid(1,1); %#ok<PFBNS>
    
    for j = 1:height(positionGrid)
        CoM{i}(2) = CoM{i}(2) + (mass{i}(j) * positionGrid(j,2));
        CoM{i}(3) = CoM{i}(3) + (mass{i}(j) * positionGrid(j,3));
    end
    
    CoM{i}(2) = CoM{i}(2) / sum(mass{i});
    CoM{i}(3) = CoM{i}(3) / sum(mass{i});
    
    send(dQ, []);
end
clear mass positionGrid;

delete(wB);

mapData.inst.CoM = CoM;
clear CoM;

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
mapData.mean.massNorm = mapData.mean.mass / massNormalisation;
clear nParticlesMean d10Mean d20Mean d30Mean d32Mean massMean;

% Calculate Time-Averaged Centre of Mass
mapData.mean.CoM = zeros(1,3);
mapData.mean.CoM(1) = mapData.positionGrid(1,1);

for i = 1:height(mapData.positionGrid)
    mapData.mean.CoM(2) = mapData.mean.CoM(2) + ...
                          (mapData.mean.mass(i) * mapData.positionGrid(i,2));
    mapData.mean.CoM(3) = mapData.mean.CoM(3) + ...
                          (mapData.mean.mass(i) * mapData.positionGrid(i,3));
end

mapData.mean.CoM(2) = mapData.mean.CoM(2) / sum(mapData.mean.mass);
mapData.mean.CoM(3) = mapData.mean.CoM(3) / sum(mapData.mean.mass);

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
    selection = input('Plot Time-Averaged Data? [y/n]: ', 's');

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
    selection = input('Plot Instantaneous Data? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotInst = false;
        
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotInst = true;
        
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

% Define Plot Limits
switch format
    
    case 'A'
        
        if contains(caseName, ["Run_Test", "Windsor"])
            xLimsPlot = [0.31875; 1.52725];
            yLimsPlot = [-0.2445; 0.2445];
            zLimsPlot = [0; 0.389];
        end
        
    case 'B'
        
        if contains(caseName, ["Run_Test", "Windsor"])
            xLimsPlot = [0.31875; 4.65925];
            yLimsPlot = [-0.5945; 0.5945];
            zLimsPlot = [0; 0.739];
        end
        
end

if normalise
    xLimsPlot = round((xLimsPlot / 1.044), spacePrecision);
    yLimsPlot = round((yLimsPlot / 1.044), spacePrecision);
    zLimsPlot = round((zLimsPlot / 1.044), spacePrecision);
end

positionData = mapData.positionGrid;
cMap = viridis(24);
figTitle = '-'; % Leave Blank ('-') for Formatting Purposes

if plotMean
    
    for i = 1:height(plotVars)
        disp(['    Presenting Time-Averaged ''', plotVars{i}, ''' Data...']);
        
        contaminantData = mapData.mean.(plotVars{i});
        
        switch format
            
            case 'A'
                figName = ['Time_Averaged_Base_', plotVars{i}, '_Map'];
                
            case 'B'
                figName = ['Time_Averaged_', planePos, '_', plotVars{i}, '_Map'];
                
        end
        
        if contains(plotVars{i}, ["mass", "massNorm"])
            CoM = mapData.mean.CoM;
        else
            CoM = [];
        end
            
        figSubtitle = ' ';
        
        if contains(plotVars{i}, ["d10", "d20", "d30", "d32"])
%             cLims = dLims;
%             cLims = [0; 30]; % Time-Averaged Base Contamination
            cLims = [0; 40]; % Time-Averaged Planar Contamination
        elseif strcmp(plotVars{i}, 'massNorm')
            cLims = [0; 1]; % Time-Averaged Base Contamination
%             cLims = [0; 20]; % Time-Averaged Planar Contaminatyion
        else
            cLims = [0; max(contaminantData)];
        end
        
        fig = contaminantPlots(xLimsPlot, yLimsPlot, zLimsPlot, xLimsData, yLimsData, zLimsData, ...
                               mapPerim, positionData, contaminantData, fig, figName, cMap, geometry, ...
                               xDims, CoM, figTitle, figSubtitle, cLims, normalise);
    end
    
    disp(' ');
end

if plotInst
    
    for i = 1:height(plotVars)
        disp(['    Presenting Instantaneous ''', plotVars{i}, ''' Data...']);
        
        figHold = fig;
        
        for j = 1:height(mapData.inst.time)
            
            if j ~= 1
                clf(fig);
                fig = figHold;
            end
            
            contaminantData = mapData.inst.(plotVars{i}){j};
            figTime = num2str(mapData.inst.time(j), ['%.', num2str(timePrecision), 'f']);
            
            switch format
                
                case 'A'
                    figName = ['Instantaneous_Base_', plotVars{i}, '_Map_T', erase(figTime, '.')];
                
                case 'B'
                    figName = ['Instantaneous_', planePos, '_', plotVars{i}, '_Map_T', erase(figTime, '.')];
            end
            
        if contains(plotVars{i}, ["mass", "massNorm"])
            CoM = mapData.mean.CoM;
        else
            CoM = [];
        end
            
        figSubtitle = ' ';
        
        if contains(plotVars{i}, ["d10", "d20", "d30", "d32"])
            cLims = dLims;
%             cLims = [0; 1]; % Time-Averaged Base Contamination
%             cLims = [0; 1]; % Time-Averaged Planar Contamination
        elseif strcmp(plotVars{i}, 'massNorm')
            cLims = [0; 1]; % Time-Averaged Base Contamination
%             cLims = [0; 1]; % Time-Averaged Planar Contaminatyion
        else
            cLims = [0; max(contaminantData)];
        end
            
            fig = contaminantPlots(xLimsPlot, yLimsPlot, zLimsPlot, xLimsData, yLimsData, zLimsData, ...
                                   mapPerim, positionData, contaminantData, fig, figName, cMap, geometry, ...
                                   xDims, CoM, figTitle, figSubtitle, cLims, normalise);
        end
        
    end
    
    disp(' ');
end

if ~plotMean && ~plotInst
    disp('    Skipping Data Presentation');
end


%% Save Map Data

valid = false;
while ~valid
    disp(' ');
    selection = input('Save Data for Future Use? [y/n]: ', 's');
    
    if selection == 'n' | selection == 'N' %#ok<OR2>
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        
        switch format
            
            case 'A'
                
                if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/contaminantMap/Base/', caseName], 'dir')
                    mkdir(['/mnt/Processing/Data/Numerical/MATLAB/contaminantMap/Base/', caseName]);
                end
                
            case 'B'
                
                if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/contaminantMap/', planePos, '/', caseName], 'dir')
                    mkdir(['/mnt/Processing/Data/Numerical/MATLAB/contaminantMap/', planePos, '/', caseName]);
                end
                
        end
        
        minD = num2str(dLims(1));
        maxD = num2str(dLims(2));
        startInst = erase(num2str(mapData.inst.time(1), ['%.', num2str(timePrecision), 'f']), '.');
        endInst = erase(num2str(mapData.inst.time(end), ['%.', num2str(timePrecision), 'f']), '.');
        
        freq = num2str(round((1 / (deltaT * sampleInterval)), timePrecision));
        
        if normalise
            fileName = ['/D', minD, '_D', maxD, '_T', startInst, '_T', endInst, '_F', freq, '_Norm.mat'];
        else
            fileName = ['/D', minD, '_D', maxD, '_T', startInst, '_T', endInst, '_F', freq, '.mat'];
        end
        
        switch format
            
            case 'A'
                save(['/mnt/Processing/Data/Numerical/MATLAB/contaminantMap/Base/', caseName, fileName], ...
                     'mapData', 'sampleInterval', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp(['    Saving to: ~/Data/Numerical/MATLAB/contaminantMap/Base/', caseName, fileName]);
                disp('        Success');
                 
            case 'B'
                save(['/mnt/Processing/Data/Numerical/MATLAB/contaminantMap/', planePos, '/', caseName, fileName], ...
                     'mapData', 'sampleInterval', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp(['    Saving to: ~/Data/Numerical/MATLAB/contaminantMap/', planePos, '/', caseName, fileName]);
                disp('        Success');
        
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end


%% Local Functions

function D = inputD(type)

    D = str2double(input(['    ', type, ' Diameter of Interest [', char(956), 'm]: '], 's'));
    
    if isnan(D) || length(D) > 1 || D < 1
        disp('        WARNING: Invalid Entry');
        D = -1;
    end
    
end
