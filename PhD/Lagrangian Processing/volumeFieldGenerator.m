%% Lagrangian Volume Field Generator v3.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

normalise = true; % Normalisation of Dimensionsn

cloudName = 'kinematicCloud'; % OpenFOAM Cloud Name

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

cellSize = 8e-3; % Spatial Resolution of Contaminant Map [m or l]

massNormalisation = 3.944150754311134e-10; % Square-Back Base Time-Average

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp ('===========================');
disp ('Volume Field Generator v3.0');
disp ('===========================');

disp (' ');
disp (' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Updated calls to 'globalPos' to 'positionCartesian'
% v1.2 - Updated to Support Changes to 'timeDirectories.m'
% v2.0 - Rewritten to Follow Recent Lagrangian Processing Structure Changes
% v2.1 - Added Time-Averaging Functionality
% v3.0 - Rewrite, Accommodating New OpenFOAM Data Formats


%% Initialise Case

[caseFolder, caseName, timeDirs, deltaT, timePrecision, geometry, ...
 xDims, yDims, zDims, spacePrecision, normalise] = initialiseCaseData(normalise);

disp(' ');
disp(' ');


%% Select Region of Interest

disp('Region of Interest');
disp('-------------------');

disp(' ');

disp('Possible Regions of Interest:');
disp('    A: Recirculation Region');
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


%% Initialise Lagrangian Data

[dataID, LagProps, ~, ~, LagData, sampleInterval] = initialiseLagData(caseFolder, caseName, ...
                                                                      cloudName, false, false, ...
                                                                      true, timeDirs, deltaT, ...
                                                                      timePrecision, nProc);
if normalise
    dataID = [dataID, '_Norm'];
end

disp(' ');
disp(' ');


%% Select Field Options

disp('Field Options');
disp('--------------');

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

if normalise
    dataID = insertBefore(dataID, '_Norm', ['_D', num2str(dLims(1)), '_D', num2str(dLims(2))]);
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
evalc('parpool(nProc);');

disp(' ');

disp('    Initialising...');

% Identify Empty Time Instances
i = 1;
while i <= height(LagData.time)
    
    if isempty(LagData.d{i})
        
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

% Specify Region Boundaries

switch format
    
    case 'A'
        
        if contains(caseName, ["Run_Test", "Windsor"])
            xLimsData = [0.31875; 1.075];
            yLimsData = [-0.3445; 0.3445];
            zLimsData = [0; 0.489];
        end
        
    case 'B'
        
        if contains(caseName, ["Run_Test", "Windsor"])
            xLimsData = [0.31875; 4.65925];
            yLimsData = [-0.5945; 0.5945];
            zLimsData = [0; 0.739];
        end
        
end
if normalise
    xLimsData = round((xLimsData / 1.044), spacePrecision);
    yLimsData = round((yLimsData / 1.044), spacePrecision);
    zLimsData = round((zLimsData / 1.044), spacePrecision);
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

switch format
    
    case {'A', 'B'}
        d = LagData.d;
        positionCartesian = LagData.positionCartesian;
        parfor i = 1:height(LagData.time)
            
            if positionCartesian{i} ~= -1
                index{i} = find(((d{i} * 1e6) >= dLims(1)) & ...
                                ((d{i} * 1e6) <= dLims(2)) & ...
                                (positionCartesian{i}(:,1) >= xLimsData (1)) & ...
                                (positionCartesian{i}(:,1) <= xLimsData (2)) & ...
                                (positionCartesian{i}(:,2) >= yLimsData (1)) & ...
                                (positionCartesian{i}(:,2) <= yLimsData (2)) & ...
                                (positionCartesian{i}(:,3) >= zLimsData (1)) & ...
                                (positionCartesian{i}(:,3) <= zLimsData (2))); %#ok<PFBNS>
            end
            
            send(dQ, []);
        end
        clear d positionCartesian;
        
end

delete(wB);

contaminantData.time = LagData.time;
LagData.time = []; % Free Memory

for i = 1:height(LagProps)
    contaminantData.(LagProps{i}) = cell(height(contaminantData.time),1);
end

for i = 1:height(contaminantData.time)
    
    for j = 1:height(LagProps)
        contaminantData.(LagProps{j}){i} = LagData.(LagProps{j}){i}(index{i},:);
        LagData.(LagProps{j}){i} = []; % Free Memory
    end
    
end

clear LagData;

disp(' ');

% Generate Instantaneous Volume Field
disp('    Generating Instantaneous Volume Field...');

switch format
    
    case {'A', 'B'}
        % Adjust Uniform Cell Size to Fit Region of Interest
        cellSizeX = (xLimsData (2) - xLimsData (1)) / round(((xLimsData (2) - xLimsData (1)) / cellSize));
        cellSizeY = (yLimsData (2) - yLimsData (1)) / round(((yLimsData (2) - yLimsData (1)) / cellSize));
        cellSizeZ = (zLimsData (2) - zLimsData (1)) / round(((zLimsData (2) - zLimsData (1)) / cellSize));
        
        cellVolume = cellSizeX * cellSizeY * cellSizeZ;
        
        [volData.x, volData.y, volData.z] = ndgrid(xLimsData(1):cellSizeX:xLimsData(2), ...
                                                   yLimsData(1):cellSizeY:yLimsData(2), ...
                                                   zLimsData(1):cellSizeZ:zLimsData(2));
end

volData.inst.time = contaminantData.time;

% Initialise Progress Bar
wB = waitbar(0, 'Assigning Particles to Mesh Nodes', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, height(volData.inst.time));

% Assign Particles to Volume Nodes
index = cell(height(volData.inst.time),1); % Array Position of Closest Mesh Node

totalParticles = cellfun(@height, contaminantData.positionCartesian);
positionCartesian = contaminantData.positionCartesian;
x = volData.x;
y = volData.y;
z = volData.z;
parfor i = 1:height(volData.inst.time)
    
    if positionCartesian{i} ~= -1
        index{i} = zeros(height(positionCartesian{i}),3);
        
        for j = 1:totalParticles(i)
            [~, index{i}(j,1)] = min(abs(positionCartesian{i}(j,1) - x(:,1,1))); %#ok<PFBNS>
            [~, index{i}(j,2)] = min(abs(positionCartesian{i}(j,2) - y(1,:,1))); %#ok<PFBNS>
            [~, index{i}(j,3)] = min(abs(positionCartesian{i}(j,3) - z(1,1,:))); %#ok<PFBNS>

        end
    
    end
    
    send(dQ, []);
end
clear totalParticles positionCartesian x y z;

delete(wB);

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Field Variables', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, height(volData.inst.time));

% Calculate Instantaneous Field Variables
nParticles = cell(height(volData.inst.time),1); % Number of Particles in Cell
d10 = nParticles; % Arithmetic Mean Diameter in Cell
d20 = nParticles; % Surface Mean Diameter in Cell
d30 = nParticles; % Volume Mean Diameter in Cell
d32 = nParticles; % Sauter Mean Diameter in Cell
mass = nParticles; % Total Mass in Cell
volFraction = nParticles; % Fraction of Cell Volume Occupied by Spray

totalParticles = cellfun(@height, contaminantData.positionCartesian);
x = volData.x;
nParticle = contaminantData.nParticle;
d = contaminantData.d;
parfor i = 1:height(volData.inst.time)
    nParticles{i} = zeros(size(x));
    d10{i} = nParticles{i};
    d20{i} = nParticles{i};
    d30{i} = nParticles{i};
    d32{i} = nParticles{i};
    mass{i} = nParticles{i};
    volFraction{i} = nParticles{i};
    
    for j = 1:totalParticles(i)
        nParticles{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = nParticles{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + ...
                                                                     nParticle{i}(j);
                                                                 
        d10{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = d10{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + ...
                                                              (nParticle{i}(j) * d{i}(j));
                                                          
        d20{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = d20{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + ...
                                                              (nParticle{i}(j) * (d{i}(j)^2));
                                                          
        d30{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = d30{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + ...
                                                              (nParticle{i}(j) * (d{i}(j)^3));
                                                          
        mass{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = mass{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + ...
                                                               (nParticle{i}(j) * ((1 / 12) * tau * (d{i}(j)^3)));
    end
    
    volFraction{i} = mass{i} / cellVolume;
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
clear totalParticles x nParticle d;

delete(wB);

clear contaminantData;

volData.inst.nParticles = nParticles;
volData.inst.d10 = d10;
volData.inst.d20 = d20;
volData.inst.d30 = d30;
volData.inst.d32 = d32;
volData.inst.mass = mass;
volData.inst.volFraction = volFraction;
clear nParticles d10 d20 d30 d32 mass volFraction;

disp(' ');

% Generate Time-Averaged Volume Field
disp('    Generating Time-Averaged Volume Field...');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Time-Averaged Field Variables', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, height(volData.inst.time));

% Calculate Time-Averaged Field Variables
nParticlesMean = zeros(size(volData.x));
d10Mean = nParticlesMean;
d20Mean = nParticlesMean;
d30Mean = nParticlesMean;
d32Mean = nParticlesMean;
massMean = nParticlesMean;
volFractionMean = nParticlesMean;

nParticles = volData.inst.nParticles;
d10 = volData.inst.d10;
d20 = volData.inst.d20;
d30 = volData.inst.d30;
d32 = volData.inst.d32;
mass = volData.inst.mass;
volFraction = volData.inst.volFraction;
parfor i = 1:height(volData.inst.time)
    nParticlesMean = nParticlesMean + nParticles{i};
    d10Mean = d10Mean + d10{i};
    d20Mean = d20Mean + d20{i};
    d30Mean = d30Mean + d30{i};
    d32Mean = d32Mean + d32{i};
    massMean = massMean + mass{i};
    volFractionMean = volFractionMean + volFraction{i};
    
    send(dQ, []);
end
clear nParticles d10 d20 d30 d32 mass volFraction;

delete(wB);

volData.mean.nParticles = nParticlesMean / height(volData.inst.time);
volData.mean.d10 = d10Mean / height(volData.inst.time);
volData.mean.d20 = d20Mean / height(volData.inst.time);
volData.mean.d30 = d30Mean / height(volData.inst.time);
volData.mean.d32 = d32Mean / height(volData.inst.time);
volData.mean.mass = massMean / height(volData.inst.time);
volData.mean.volFraction = volFractionMean / height(volData.inst.time);
clear nParticlesMean d10Mean d20Mean d30Mean d32Mean massMean volFractionMean;

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


%% Local Functions

function D = inputD(type)

    D = str2double(input(['    ', type, ' Diameter of Interest [', char(956), 'm]: '], 's'));
    
    if isnan(D) || length(D) > 1 || D < 1
        disp('        WARNING: Invalid Entry');
        D = -1;
    end
    
end