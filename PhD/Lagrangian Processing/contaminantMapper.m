%% Lagrangian Contamination Mapper v2.0

clear variables;
close all;
clc;

normalise = true; % Normalisation of Dimensions

cloudName = 'kinematicCloud'; % OpenFOAM Cloud Name

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

cellSize = 8e-3; % Spatial Resolution of Contaminant Map [m or l]

massNormalisation = 3.744918231958561e-10; % Square-Back Base Time-Average

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
        [LagProps, ~, LagData, ~] = initialiseLagData(caseFolder, caseName, cloudName, ...
                                                      false, true, false, ...
                                                      timeDirs, deltaT, timePrecision, nProc);
                                                                                
    case 'B'
        [LagProps, LagData, ~, ~] = initialiseLagData(caseFolder, caseName, cloudName, ...
                                                      true, false, false, ...
                                                      timeDirs, deltaT, timePrecision, nProc);

        % Select Plane of Interest
        planes = fieldnames(LagData);
        
        valid = false;
        while ~valid
            [index, valid] = listdlg('listSize', [300, 300], ...
                                     'selectionMode', 'single', ...
                                     'name', 'Select Plane of Interest', ...
                                     'listString', planes);
            
            if ~valid
                disp('    WARNING: No Case Type Selected');
            end
        
        end
        clear valid;

        LagData = LagData.(planes{index});
        clear planes;
end

disp(' ');
disp(' ');


%% Select Mapping Options

disp('Mapping Options');
disp('----------------');

disp(' ');

disp('Loaded Contaminant Data Spans:');

disp(['    T = ', num2str(LagData.time(1), ['%.', num2str(timePrecision), 'f']), ' s ', ...
       '-> T = ', num2str(LagData.time(end), ['%.', num2str(timePrecision), 'f']), ' s']);
disp(['        ', char(916), 'T = ' num2str(deltaT, ['%.', num2str(timePrecision), 'f']), 's']);

valid = false;
while ~valid
    disp(' ');
    selection = input('Utilise All Available Time Instances? [y/n]: ', 's');
    
    if selection == 'n' | selection == 'N' %#ok<OR2>
        timeInsts = inputTimes(LagData.time);
        
        if timeInsts == -1
            continue
        elseif isempty(timeInsts)
            disp('        WARNING: No Lagrangian Data Available for Selected Time Instances');
            continue
        end
        
        % Remove Unnecessary Data
        LagData.time = LagData.time(timeInsts);
        LagData.timeExact = LagData.timeExact(timeInsts);
        
        for i = 1:height(LagProps)
            LagData.(LagProps{i}) = LagData.(LagProps{i})(timeInsts);
        end
        
        valid = true;        
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

% Remove Empty Time Instances
i = 1;
while i <= height(LagData.time)
    
    if isempty(LagData.timeExact{i})
        LagData.time(i) =[];
        LagData.timeExact(i) = [];
        
        for j = 1:height(LagProps)
            LagData.(LagProps{j})(i) = [];
        end
        
    else
        i = i + 1;
    end
    
end
clear i;

dLims = zeros(2,1);

valid = false;
while ~valid
    disp(' ');
    selection = input('Filter Particle Diameters? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        dLims = [0; 120];
        
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
        
        if (dLims(2) < 0) || (dLims(1) > 120)
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

% Shift Data Origin
if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
    
    for i = 1:height(LagData.time)
        LagData.positionCartesian{i}(:,1) = LagData.positionCartesian{i}(:,1) + 1.325;
    end
    
end

% Normalise Dimensions
if normalise
    
    if contains(caseName, ["Run_Test", "Windsor"])
        
        for i = 1:height(LagData.time)
            LagData.positionCartesian{i}  = round((LagData.positionCartesian{i} / 1.044), spacePrecision);
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
        
        basePerim = boundary(basePoints(:,2), basePoints(:,3), 0.95);
        basePerim = basePoints(basePerim,:);
        basePoly = polyshape(basePerim(:,2), basePerim(:,3), 'keepCollinearPoints', true);
        basePoly = polybuffer(basePoly, -0.0025, 'jointType', 'square');
        basePerim = ones(height(basePoly.Vertices),3) * basePerim(1,1);
        basePerim(:,[2,3]) = basePoly.Vertices(:,[1,2]);

        if ~all(basePerim(1,:) == basePerim(end,:))
            basePerim = vertcat(basePerim, basePerim(1,:)); % Close Boundary
        end
        
        clear basePoly;
        
        xLims = xDims(2);
        yLims = [min(basePoints(:,2)); max(basePoints(:,2))];
        zLims = [min(basePoints(:,3)); max(basePoints(:,3))];
        
    case 'B'
        if contains(caseName, 'Run_Test') || contains(caseName, 'Windsor')
            basePerim = [];
            
            xLims = LagData.positionCartesian{end}(1,1);
            yLims = [-0.5945; 0.5945];
            zLims = [0; 0.739];
            
            if normalise
                yLims = round((yLims / 1.044), spacePrecision);
                zLims = round((zLims / 1.044), spacePrecision);
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
    index{i} = find(((d{i} * 1e6) >= dLims(1)) & ...
                    ((d{i} * 1e6) <= dLims(2)) & ...
                    (positionCartesian{i}(:,1) == xLims) & ...
                    (positionCartesian{i}(:,2) >= yLims(1)) & ...
                    (positionCartesian{i}(:,2) <= yLims(2)) & ...
                    (positionCartesian{i}(:,3) >= zLims(1)) & ...
                    (positionCartesian{i}(:,3) <= zLims(2))); %#ok<PFBNS>
    
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
disp('    Generating Instantaneous Contaminant Maps');

cellSizeX = cellSize;
cellSizeY = (yLims(2) - yLims(1)) / round(((yLims(2) - yLims(1)) / cellSize));
cellSizeZ = (zLims(2) - zLims(1)) / round(((zLims(2) - zLims(1)) / cellSize));

[y, z] = meshgrid(yLims(1):cellSizeY:yLims(2), zLims(1):cellSizeZ:zLims(2));

mapData.positionGrid = zeros(height(y(:)),3);
mapData.positionGrid(:,1) = xLims;
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
    index{i} = dsearchn(positionGrid, positionCartesian{i});
    
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

positionGrid = mapData.positionGrid;
positionCartesian = contaminantData.positionCartesian;
nParticle = contaminantData.nParticle;
d = contaminantData.d;
d32_tmp = nParticles;
parfor i = 1:height(mapData.inst.time)
    nParticles{i} = zeros(height(positionGrid),1);
    d10{i} = nParticles{i};
    d20{i} = nParticles{i};
    d30{i} = nParticles{i};
    d32{i} = nParticles{i};
    mass{i} = nParticles{i};
    
    d32_tmp{i} = nParticles{i};
    
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
clear nParticles d10 d20 d30 d32 mass;

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

disp(['    Read Time: ', num2str(executionTime), 's']);

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
    % Select Variables of Interest
    plotVars = fieldnames(mapData.mean);
    plotVars = plotVars(1:(end - 1));

    valid = false;
    while ~valid
        [index, valid] = listdlg('listSize', [300, 300], ...
                                 'selectionMode', 'multiple', ...
                                 'name', 'Select Variable(s) to Plot', ...
                                 'listString', plotVars);

        if ~valid
            disp('    WARNING: No Mapping Variables Selected');
        end
    end
    clear valid;

    plotVars = plotVars(index);
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
            yLimsPlot = yLims;
            zLimsPlot = zLims;
        end
        
end

if normalise
    xLimsPlot = round((xLimsPlot / 1.044), spacePrecision);
    yLimsPlot = round((yLimsPlot / 1.044), spacePrecision);
    zLimsPlot = round((zLimsPlot / 1.044), spacePrecision);
end

xLimsData = xLims;
yLimsData = yLims;
zLimsData = zLims;
positionData = mapData.positionGrid;
cMap = viridis(24);
figTitle = '-'; % Leave Blank ('-') for Formatting Purposes

if plotMean
    
    for i = 1:height(plotVars)
        disp(['    Presenting Time-Averaged ''', plotVars{i}, ''' Data...']);
        
        contaminantData = mapData.mean.(plotVars{i});
        
        switch format
            
            case 'A'
                figName = ['Time_Averaged_Base', plotVars{i}, '_Map'];
                
            case 'B'
                figName = ['Time_Averaged_Planar', plotVars{i}, '_Map'];
                
        end
        
        CoM = mapData.mean.CoM;
        figSubtitle = ' ';
        
        if contains(plotVars{i}, ["d10", "d20", "d30", "d32"])
            cLims = dLims;
        elseif strcmp(plotVars{i}, 'massNorm')
            cLims = [0; 1];
        else
            cLims = [0; max(contaminantData)];
        end
        
        fig = contaminantPlots(xLimsPlot, yLimsPlot, zLimsPlot, xLimsData, yLimsData, zLimsData, ...
                               basePerim, positionData, contaminantData, fig, figName, cMap, geometry, ...
                               xDims, CoM, figTitle, figSubtitle, cLims, normalise);
    end
    
    disp(' ');
end

if plotInst
    
    for i = 1:height(plotVars)
        disp(['    Presenting Instantaneous ''', plotVars{i}, ''' Data...']);
        
        for j = 1:height(mapData.inst.time)
            contaminantData = mapData.inst.(plotVars{i}){j};
            figTime = num2str(mapData.inst.time{j}, ['%.', num2str(timePrecision), 'f']);
            
            switch format
                
                case 'A'
                    figName = ['Instantaneous_Base', plotVars{i}, '_Map_T', figTime];
                
                case 'B'
                    figName = ['Instantaneous_Planar', plotVars{i}, '_Map_T', figTime];
            end
            
            CoM = mapData.inst.CoM{j};
            figSubtitle = figTime;
            
            if contains(plotVars{i}, ["d10", "d20", "d30", "d32"])
                cLims = dLims;
            elseif strcmp(plotVars{i}, 'massNorm')
                cLims = [0; 1];
            else
                cLims = [0; max(contaminantData)];
            end
            
            fig = contaminantPlots(xLimsPlot, yLimsPlot, zLimsPlot, xLimsData, yLimsData, zLimsData, ...
                                   basePerim, positionData, contaminantData, fig, figName, cMap, geometry, ...
                                   xDims, CoM, figTitle, figSubtitle, cLims, normalise);
        end
        
    end
    
    disp(' ');
end

if ~plotMean && ~plotInst
    disp('    Skipping Data Presentation');
end


%% Save Map Data

% Save Shit


%% Local Functions

function timeInsts = inputTimes(origTimes)

    timeInsts = str2num(input('    Input Desired Time Instances (Row Vector Form) [s]: ', 's')); %#ok<ST2NM>
    
    if any(isnan(timeInsts)) || ~isrow(timeInsts) > 1 || any(timeInsts <= 0)
        disp('        WARNING: Invalid Entry');
        timeInsts = -1;
    else
        timeInsts = find(ismember(origTimes, timeInsts));
    end

end


function D = inputD(type)

    D = str2double(input(['    ', type, ' Diameter of Interest [', char(956), 'm]: '], 's'));
    
    if isnan(D) || length(D) > 1
        disp('        WARNING: Invalid Entry');
        D = -1;
    end
    
end