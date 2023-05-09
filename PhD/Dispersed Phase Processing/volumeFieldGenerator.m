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


%% Lagrangian Volume Field Generator v3.2

cellSize.target = 8e-3; % Target Spatial Resolution of Volume Field [m]

cloudName = 'kinematicCloud'; % OpenFOAM Cloud Name

normalise = false; % Normalisation of Dimensions


disp('===========================');
disp('Volume Field Generator v3.2');
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
% v3.1 - Added Support for Multiple Mean Particle Diameter Definitions
% v3.2 - Removed Parallelisation of Field Calculations to Reduce Memory Requirements


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


%% Initialise Lagrangian Data

[dataID, LagProps, ~, ~, LagData, sampleInterval] = initialiseLagData(saveLocation, caseFolder, caseName, ...
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
            LagData.(LagProps{j}){i} = [];
        end
        
    else
        i = i + 1;
    end
    
end

nTimes = height(LagData.time);

% Shift Data Origin
if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
    
    for i = 1:nTimes
        
        if ~isempty(LagData.positionCartesian{i})
            LagData.positionCartesian{i}(:,1) = LagData.positionCartesian{i}(:,1) + 1.325;
        end
        
    end
    
end

% Normalise Dimensions
if normalise
    
    if contains(caseName, ["Run_Test", "Windsor"])
        
        for i = 1:nTimes
            
            if ~isempty(LagData.positionCartesian{i})
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
            xLimsData = [0.31875; 1.26625];
            yLimsData = [-0.3445; 0.3445];
            zLimsData = [0; 0.489];
        end
        
    case 'B'
        
        if contains(caseName, ["Run_Test", "Windsor"])
            xLimsData = [0.31875; 2.57125]; % 2L
            yLimsData = [-0.5945; 0.5945];
            zLimsData = [0; 0.639];
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

parforWaitBar(wB, nTimes);

% Collate Particles of Interest
index = cell(height(LagData.time),1);

d = LagData.d;
positionCartesian = LagData.positionCartesian;
parfor i = 1:height(LagData.time)
    
    if ~isempty(positionCartesian{i})
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
        
delete(wB);

% Remove Unnecessary Data
LagData = rmfield(LagData, {'origId', 'origProcId', 'U'});
LagProps = {'d'; 'nParticle'; 'positionCartesian'};

for i = 1:nTimes
    
    for j = 1:height(LagProps)
        LagData.(LagProps{j}){i} = LagData.(LagProps{j}){i}(index{i},:);
    end
    
end

disp(' ');

% Generate Instantaneous Volume Field
disp('    Generating Instantaneous Volume Field...');

if normalise
    cellSize.target = round((cellSize.target / 1.044), spacePrecision);
end

% Adjust Uniform Cell Size to Fit Region of Interest
cellSize.x = (xLimsData(2) - xLimsData(1)) / round(((xLimsData(2) - xLimsData(1)) / cellSize.target));
cellSize.y = (yLimsData(2) - yLimsData(1)) / round(((yLimsData(2) - yLimsData(1)) / cellSize.target));
cellSize.z = (zLimsData(2) - zLimsData(1)) / round(((zLimsData(2) - zLimsData(1)) / cellSize.target));

cellSize.volume = cellSize.x * cellSize.y * cellSize.z;

[volumeData.x, volumeData.y, volumeData.z] = ndgrid(xLimsData(1):(cellSize.x):xLimsData(2), ...
                                                    yLimsData(1):(cellSize.y):yLimsData(2), ...
                                                    zLimsData(1):(cellSize.z):zLimsData(2));

volumeData.inst.time = LagData.time;

% Initialise Progress Bar
wB = waitbar(0, 'Assigning Particles to Mesh Nodes', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, nTimes);

% Assign Particles to Volume Nodes
totalParticles = cellfun(@height, LagData.positionCartesian);
index = cell(nTimes,1); % Array Position of Closest Mesh Node

positionCartesian = LagData.positionCartesian;
x = volumeData.x;
y = volumeData.y;
z = volumeData.z;
parfor i = 1:nTimes
    
    if totalParticles(i) > 0
        index{i} = zeros(height(positionCartesian{i}),3);
        
        for j = 1:totalParticles(i)
            [~, index{i}(j,1)] = min(abs(positionCartesian{i}(j,1) - x(:,1,1))); %#ok<PFBNS>
            [~, index{i}(j,2)] = min(abs(positionCartesian{i}(j,2) - y(1,:,1))); %#ok<PFBNS>
            [~, index{i}(j,3)] = min(abs(positionCartesian{i}(j,3) - z(1,1,:))); %#ok<PFBNS>

        end
    
    end
    
    send(dQ, []);
end
clear positionCartesian x y z;

delete(wB);

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Field Variables', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Calculate Instantaneous Field Variables
volumeData.inst.nParticles = cell(nTimes,1);    % Number of Particles in Cell
volumeData.inst.volFraction = cell(nTimes,1);   % Fraction of Cell Volume Occupied by Spray
volumeData.inst.mass = cell(nTimes,1);          % Total Mass in Cell
% volumeData.inst.d43 = cell(nTimes,1);           % De Brouckere Mean Diameter in Cell
% volumeData.inst.d32 = cell(nTimes,1);           % Sauter Mean Diameter in Cell
% volumeData.inst.d30 = cell(nTimes,1);           % Volume Mean Diameter in Cell
% volumeData.inst.d20 = cell(nTimes,1);           % Surface Mean Diameter in Cell
volumeData.inst.d10 = cell(nTimes,1);           % Arithmetic Mean Diameter in Cell

totalParticles = cellfun(@height, LagData.positionCartesian);
for i = 1:nTimes
    volumeData.inst.nParticles{i} = zeros(size(volumeData.x));
    volumeData.inst.volFraction{i} = volumeData.inst.nParticles{i};
    volumeData.inst.mass{i} = volumeData.inst.nParticles{i};
%     volumeData.inst.d43{i} = volumeData.inst.nParticles{i};
%     volumeData.inst.d32{i} = volumeData.inst.nParticles{i};
%     volumeData.inst.d30{i} = volumeData.inst.nParticles{i};
%     volumeData.inst.d20{i} = volumeData.inst.nParticles{i};
    volumeData.inst.d10{i} = volumeData.inst.nParticles{i};
    
    if totalParticles(i) > 0

        for j = 1:totalParticles(i)
            volumeData.inst.nParticles{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = ...
            volumeData.inst.nParticles{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + LagData.nParticle{i}(j);
    
            volumeData.inst.mass{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = ...
            volumeData.inst.mass{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + (LagData.nParticle{i}(j) * ((1 / 12) * tau * (LagData.d{i}(j)^3)));
    
%             volumeData.inst.d43{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = ...
%             volumeData.inst.d43{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + (LagData.nParticle{i}(j) * (LagData.d{i}(j)^4));
    
%             volumeData.inst.d30{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = ...
%             volumeData.inst.d30{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + (LagData.nParticle{i}(j) * (LagData.d{i}(j)^3));
    
%             volumeData.inst.d20{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = ...
%             volumeData.inst.d20{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + (LagData.nParticle{i}(j) * (LagData.d{i}(j)^2));
                                                                     
            volumeData.inst.d10{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = ...
            volumeData.inst.d10{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + (LagData.nParticle{i}(j) * LagData.d{i}(j));
        end

    end

    % Remove Unnecessary Data
    index{i} = [];
    LagData.nParticle{i} = [];
    LagData.d{i} = [];
    
    % Calculate Derived Variables
    volumeData.inst.volFraction{i} = volumeData.inst.mass{i} / cellSize.volume;
    volumeData.inst.mass{i} = 1000 * volumeData.inst.mass{i};
%     volumeData.inst.d43{i} = (volumeData.inst.d43{i} ./ volumeData.inst.d30{i}) * 1e6;
%     volumeData.inst.d32{i} = (volumeData.inst.d30{i} ./ volumeData.inst.d20{i}) * 1e6;
%     volumeData.inst.d30{i} = ((volumeData.inst.d30{i} ./ volumeData.inst.nParticles{i}).^(1/3)) * 1e6;
%     volumeData.inst.d20{i} = ((volumeData.inst.d20{i} ./ volumeData.inst.nParticles{i}).^(1/2)) * 1e6;
    volumeData.inst.d10{i} = (volumeData.inst.d10{i} ./ volumeData.inst.nParticles{i}) * 1e6;
    
    % Set Empty Cells Back to Zero
    volumeData.inst.nParticles{i}(isnan(volumeData.inst.nParticles{i})) = 0;
    volumeData.inst.volFraction{i}(isnan(volumeData.inst.volFraction{i})) = 0;
    volumeData.inst.mass{i}(isnan(volumeData.inst.mass{i})) = 0;
%     volumeData.inst.d43{i}(isnan(volumeData.inst.d43{i})) = 0;
%     volumeData.inst.d32{i}(isnan(volumeData.inst.d32{i})) = 0;
%     volumeData.inst.d30{i}(isnan(volumeData.inst.d30{i})) = 0;
%     volumeData.inst.d20{i}(isnan(volumeData.inst.d20{i})) = 0;
    volumeData.inst.d10{i}(isnan(volumeData.inst.d10{i})) = 0;
    
    waitbar((i / nTimes), wB);
end

delete(wB);

clear LagData;

disp(' ');

% Generate Time-Averaged Volume Field
disp('    Generating Time-Averaged Volume Field...');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Time-Averaged Field Variables', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Calculate Time-Averaged Field Variables
volumeData.mean.nParticles = zeros(size(volumeData.x));
volumeData.mean.volFraction = volumeData.mean.nParticles;
volumeData.mean.mass = volumeData.mean.nParticles;
% volumeData.mean.d43 = volumeData.mean.nParticles;
% volumeData.mean.d32 = volumeData.mean.nParticles;
% volumeData.mean.d30 = volumeData.mean.nParticles;
% volumeData.mean.d20 = volumeData.mean.nParticles;
volumeData.mean.d10 = volumeData.mean.nParticles;

for i = 1:nTimes
    volumeData.mean.nParticles = volumeData.mean.nParticles + volumeData.inst.nParticles{i};
    volumeData.mean.volFraction = volumeData.mean.volFraction + volumeData.inst.volFraction{i};
    volumeData.mean.mass = volumeData.mean.mass + volumeData.inst.mass{i};
%     volumeData.mean.d43 = volumeData.mean.d43 + volumeData.inst.d43{i};
%     volumeData.mean.d32 = volumeData.mean.d32 + volumeData.inst.d32{i};
%     volumeData.mean.d30 = volumeData.mean.d30 + volumeData.inst.d30{i};
%     volumeData.mean.d20 = volumeData.mean.d20 + volumeData.inst.d20{i};
    volumeData.mean.d10 = volumeData.mean.d10 + volumeData.inst.d10{i};
    
    waitbar((i / nTimes), wB);
end

volumeData.mean.nParticles = volumeData.mean.nParticles / nTimes;
volumeData.mean.volFraction = volumeData.mean.volFraction / nTimes;
volumeData.mean.mass = volumeData.mean.mass / nTimes;
% volumeData.mean.d43 = volumeData.mean.d43 / nTimes;
% volumeData.mean.d32 = volumeData.mean.d32 / nTimes;
% volumeData.mean.d30 = volumeData.mean.d30 / nTimes;
% volumeData.mean.d20 = volumeData.mean.d20 / nTimes;
volumeData.mean.d10 = volumeData.mean.d10 / nTimes;

delete(wB);

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
        nFrames = inputFrames(nTimes);
        
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


%% Present Volume Fields

disp('Volume Field Presentation');
disp('--------------------------');

disp(' ');

if plotInst || plotMean
    xInit = volumeData.x;
    yInit = volumeData.y;
    zInit = volumeData.z;
    POD = false;
    
    if contains(caseName, 'Windsor')
        cMap = viridis(3);
        
        if strcmp(caseName, 'Windsor_SB_wW_Upstream_SC')
            cMap = cMap(1,:);
        elseif strcmp(caseName, 'Windsor_ST_20D_wW_Upstream_SC')
            cMap = cMap(2,:);
        elseif strcmp(caseName, 'Windsor_RSST_16D_U50_wW_Upstream_SC')
            cMap = cMap(3,:);
        end
        
    else
        cMap = viridis(1);
    end

    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    xLimsPlot = xLimsData;
    yLimsPlot = yLimsData;
    zLimsPlot = zLimsData;
end

if plotMean
    disp('    Presenting Time-Averaged Volume Field...');
    
    fieldData = volumeData.mean.volFraction;
%     isoValue = 1e-6; % Inner Bounds
    isoValue = 1e-8; % Outer Bounds
    
    switch format

        case 'A'
            figName = ['Near_Field_Time_Averaged_Spray_Volume_Fraction_', num2str(isoValue), ...
                       '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];

        case 'B'
            figName = ['Far_Field_Time_Averaged_Spray_Volume_Fraction_', num2str(isoValue), ...
                       '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];

    end
    
    figSubtitle = ' ';
    
    fig = plotVolumeField(xLimsData, yLimsData, zLimsData, xInit, yInit, zInit, POD, fieldData, ...
                          fig, figName, geometry, isoValue, cMap, figTitle, figSubtitle, ...
                          xLimsPlot, yLimsPlot, zLimsPlot);
                       
    disp(' ');
end

if plotInst
    disp('    Presenting Instantaneous Volume Field...');
    
    figHold = fig;
    
    isoValue = 1e-6;
    
    for i = 1:nFrames
        
        if i ~= 1
            clf(fig);
            fig = figHold;
        end
        
        fieldData = volumeData.inst.volFraction{i};
        figTime = num2str(volumeData.inst.time(i), ['%.', num2str(timePrecision), 'f']);

        switch format

            case 'A'
                figName = ['Near_Field_Instantaneous_Spray_Volume_Fraction_', num2str(isoValue), ...
                           '_T', erase(figTime, '.'), '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];

            case 'B'
                figName = ['Far_Field_Instantaneous_Spray_Volume_Fraction_', num2str(isoValue), ...
                           '_T', erase(figTime, '.'), '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];

        end
        
        figSubtitle = [figTime, ' \it{s}'];
        
        fig = plotVolumeField(xLimsData, yLimsData, zLimsData, xInit, yInit, zInit, POD, fieldData, ...
                              fig, figName, geometry, isoValue, cMap, figTitle, figSubtitle, ...
                              xLimsPlot, yLimsPlot, zLimsPlot);
                           
    end
    
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
        
        % Format volumeData to Be Human-Readable
        volumeData.positionGrid = [volumeData.x(:), volumeData.y(:), volumeData.z(:)];
        volumeData = rmfield(volumeData, {'x', 'y', 'z'});
        
        for i = 1:nTimes
            volumeData.inst.nParticles{i} = volumeData.inst.nParticles{i}(:);
            volumeData.inst.volFraction{i} = volumeData.inst.volFraction{i}(:);
            volumeData.inst.mass{i} = volumeData.inst.mass{i}(:);
%             volumeData.inst.d43{i} = volumeData.inst.d43{i}(:);
%             volumeData.inst.d32{i} = volumeData.inst.d32{i}(:);
%             volumeData.inst.d30{i} = volumeData.inst.d30{i}(:);
%             volumeData.inst.d20{i} = volumeData.inst.d20{i}(:);
            volumeData.inst.d10{i} = volumeData.inst.d10{i}(:);
        end

        volumeData.mean.nParticles = volumeData.mean.nParticles(:);
        volumeData.mean.volFraction = volumeData.mean.volFraction(:);
        volumeData.mean.mass = volumeData.mean.mass(:);
%         volumeData.mean.d43 = volumeData.mean.d43(:);
%         volumeData.mean.d32 = volumeData.mean.d32(:);
%         volumeData.mean.d30 = volumeData.mean.d30(:);
%         volumeData.mean.d20 = volumeData.mean.d20(:);
        volumeData.mean.d10 = volumeData.mean.d10(:);
        
        % Save Data
        switch format
            
            case 'A'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/volumeField/', caseName, '/nearField'], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/volumeField/', caseName, '/nearField']);
                end
                
            case 'B'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/volumeField/', caseName, '/farField'], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/volumeField/', caseName, '/farField']);
                end
                
        end
        
        switch format
            
            case 'A'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/volumeField/', caseName, '/nearField/', dataID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/volumeField/', caseName, '/nearField/', dataID, '.mat'], ...
                     'dataID', 'volumeData', 'cellSize', 'sampleInterval', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/volumeField/', caseName, '/farField/', dataID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/volumeField/', caseName, '/farField/', dataID, '.mat'], ...
                     'dataID', 'volumeData', 'cellSize', 'sampleInterval', 'dLims', 'normalise', '-v7.3', '-noCompression');
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


function nFrames = inputFrames(Nt)

    nFrames = str2double(input(['    Input Desired Frame Count [1-', num2str(Nt), ']: '], 's'));
    
    if isnan(nFrames) || nFrames <= 0 || nFrames > Nt
        disp('        WARNING: Invalid Entry');
        nFrames = -1;
    end

end
