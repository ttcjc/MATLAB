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


%% Lagrangian Volume Field Generator v3.2

cellSize.target = 8e-3; % Target Spatial Resolution of Volume Field [m]

cloudName = 'kinematicCloud'; % OpenFOAM Cloud Name

figSave = false; % Save .fig File(s)

normalise = true; % Normalisation of Dimensions

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
% v3.1 - Added Support for Multiple Mean Particle Diameter Definitions (Disabled by Default)
% v3.2 - Shifted to Thread-Based Parallelization to Reduce Memory Requirements


%% Initialise Case

[caseFolder, caseID, timeDirs, deltaT, timePrecision, geometry, ...
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
clear valid;

disp(' ');
disp(' ');


%% Initialise Lagrangian Data

[dataID, LagProps, ~, ~, LagData, sampleInterval] = initialiseLagData(saveLocation, caseFolder, caseID, ...
                                                                      cloudName, false, false, ...
                                                                      true, timeDirs, deltaT, ...
                                                                      timePrecision, nProc);

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

% Specify Region Boundaries
switch format
    
    case 'A'
        
        if contains(caseID, ["Run_Test", "Windsor"])
            xLimsData = [0.31875; 1.26625]; % 0.75L
            yLimsData = [-0.3915; 0.3915];
            zLimsData = [0; 0.4698];
        end
        
    case 'B'
        
        if contains(caseID, ["Run_Test", "Windsor"])
            xLimsData = [0.31875; 2.57125]; % 2L
            yLimsData = [-0.522; 0.522];
            zLimsData = [0; 0.6264];
        end
        
end

if normalise
    xLimsData = round((xLimsData / 1.044), spacePrecision);
    yLimsData = round((yLimsData / 1.044), spacePrecision);
    zLimsData = round((zLimsData / 1.044), spacePrecision);
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
clear i d positionCartesian;
        
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
disp('    Generating 3D Presentation Grid...');

% Adjust Uniform Cell Size to Fit Region of Interest
if normalise
    cellSize.target = round((cellSize.target / 1.044), spacePrecision);
end

cellSize.x = (xLimsData(2) - xLimsData(1)) / round(((xLimsData(2) - xLimsData(1)) / cellSize.target));
cellSize.y = (yLimsData(2) - yLimsData(1)) / round(((yLimsData(2) - yLimsData(1)) / cellSize.target));
cellSize.z = (zLimsData(2) - zLimsData(1)) / round(((zLimsData(2) - zLimsData(1)) / cellSize.target));

cellSize.volume = cellSize.x * cellSize.y * cellSize.z;

[x, y, z] = ndgrid(xLimsData(1):(cellSize.x):xLimsData(2), ...
                   yLimsData(1):(cellSize.y):yLimsData(2), ...
                   zLimsData(1):(cellSize.z):zLimsData(2));

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
index = cell(nTimes,1); % Array Position of Closest Mesh Node

positionCartesian = LagData.positionCartesian;
parfor i = 1:nTimes
    
    if totalParticles(i) > 0
        index{i} = zeros(totalParticles(i),3);
        
        for j = 1:totalParticles(i)
            [~, index{i}(j,1)] = min(abs(positionCartesian{i}(j,1) - x(:,1,1))); %#ok<PFBNS>
            [~, index{i}(j,2)] = min(abs(positionCartesian{i}(j,2) - y(1,:,1))); %#ok<PFBNS>
            [~, index{i}(j,3)] = min(abs(positionCartesian{i}(j,3) - z(1,1,:))); %#ok<PFBNS>
        end
    
    end
    
    % Remove Unnecessary Data
    positionCartesian{i} = [];
    
    % Update Waitbar
    send(dQ, []);
end
clear i positionCartesian;

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
nParticles = cell(nTimes,1); nParticles(:) = {zeros(size(x), 'single')}; % Number of Particles in Cell
volFraction = nParticles; % Fraction of Cell Volume Occupied by Spray
mass = nParticles; % Total Mass in Cell
% d32 = nParticles; % Sauter Mean Diameter in Cell
% d30 = nParticles; % Volume Mean Diameter in Cell
d20 = nParticles; % Surface Mean Diameter in Cell
d10 = nParticles; % Arithmetic Mean Diameter in Cell

nParticle = LagData.nParticle;
d = LagData.d;
cellVolume = cellSize.volume;
parfor i = 1:nTimes
    
    if totalParticles(i) > 0

        for j = 1:totalParticles(i)
            nParticles{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = ...
            nParticles{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + ...
            nParticle{i}(j);
            
            mass{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = ...
            mass{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + ...
            (nParticle{i}(j) * ((1 / 12) * tau * (d{i}(j)^3)));
            
%             d30{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = ...
%             d30{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + ...
%             (nParticle{i}(j) * (d{i}(j)^3));
            
            d20{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = ...
            d20{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + ...
            (nParticle{i}(j) * (d{i}(j)^2));
            
            d10{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = ...
            d10{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + ...
            (nParticle{i}(j) * d{i}(j));
        end

    end

    % Remove Unnecessary Data
    index{i} = [];
    nParticle{i} = [];
    d{i} = [];

    % Calculate Derived Variables
    volFraction{i} = mass{i} / cellVolume;
    mass{i} = 1000 * mass{i};
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
clear i j nParticle d cellVolume;

delete(wB);

clear LagData;

volumeData.inst.nParticles = nParticles; clear nParticles;
volumeData.inst.volFraction = volFraction; clear volFraction;
volumeData.inst.mass = mass; clear mass;
% volumeData.inst.d32 = d32; clear d32;
% volumeData.inst.d30 = d30; clear d30;
volumeData.inst.d20 = d20; clear d20;
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
volumeData.mean.nParticles = zeros(size(x), 'single');
volumeData.mean.volFraction = volumeData.mean.nParticles;
volumeData.mean.mass = volumeData.mean.nParticles;
% volumeData.mean.d32 = volumeData.mean.nParticles;
% volumeData.mean.d30 = volumeData.mean.nParticles;
volumeData.mean.d20 = volumeData.mean.nParticles;
volumeData.mean.d10 = volumeData.mean.nParticles;

for i = 1:nTimes
    volumeData.mean.nParticles = volumeData.mean.nParticles + volumeData.inst.nParticles{i};
    volumeData.mean.volFraction = volumeData.mean.volFraction + volumeData.inst.volFraction{i};
    volumeData.mean.mass = volumeData.mean.mass + volumeData.inst.mass{i};
%     volumeData.mean.d32 = volumeData.mean.d32 + volumeData.inst.d32{i};
%     volumeData.mean.d30 = volumeData.mean.d30 + volumeData.inst.d30{i};
    volumeData.mean.d20 = volumeData.mean.d20 + volumeData.inst.d20{i};
    volumeData.mean.d10 = volumeData.mean.d10 + volumeData.inst.d10{i};

    % Update Waitbar
    waitbar((i / nTimes), wB);
end
clear i;

volumeData.mean.nParticles = volumeData.mean.nParticles / nTimes;
volumeData.mean.volFraction = volumeData.mean.volFraction / nTimes;
volumeData.mean.mass = volumeData.mean.mass / nTimes;
% volumeData.mean.d32 = volumeData.mean.d32 / nTimes;
% volumeData.mean.d30 = volumeData.mean.d30 / nTimes;
volumeData.mean.d20 = volumeData.mean.d20 / nTimes;
volumeData.mean.d10 = volumeData.mean.d10 / nTimes;

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
    xInit = x;
    yInit = y;
    zInit = z;
    POD = false;
    
    if contains(caseID, 'Windsor')
        cMap = viridis(3);
        
        if strcmp(caseID, 'Windsor_SB_wW_Upstream_SC')
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

        case 'A'

            if contains(caseID, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 1.43075]; % 0.75L
                yLimsPlot = [-0.3915; 0.3915];
                zLimsPlot = [0; 0.4698];
            end

        case 'B'

            if contains(caseID, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 2.73575]; % 2L
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
    
end

if plotMean
    disp('    Presenting Time-Averaged Volume Field...');
    
    fieldData = volumeData.mean.volFraction;
    isoValue = [2e-6; 2e-8];
    figSubtitle = ' ';
    
    for i = 1:height(isoValue)

        switch format

            case 'A'
                figName = ['Near_Field_Time_Averaged_Spray_Volume_Fraction_', num2str(isoValue(i))];

            case 'B'
                figName = ['Far_Field_Time_Averaged_Spray_Volume_Fraction_', num2str(isoValue(i))];

        end

        fig = plotVolumeField(xLimsData, yLimsData, zLimsData, xInit, yInit, zInit, POD, fieldData, ...
                              fig, figName, geometry, isoValue(i), cMap, figTitle, figSubtitle, ...
                              xLimsPlot, yLimsPlot, zLimsPlot, figSave);
    end
    clear i;
                       
    disp(' ');
end

if plotInst
    disp('    Presenting Instantaneous Volume Field...');
    
    figHold = fig;
    
    isoValue = 2e-6;
    
    for i = 1:nFrames
        
        if i ~= 1
            clf(fig);
            fig = figHold;
        end
        
        fieldData = volumeData.inst.volFraction{i};
        figTime = num2str(volumeData.time(i), ['%.', num2str(timePrecision), 'f']);

        switch format

            case 'A'
                figName = ['Near_Field_Instantaneous_Spray_Volume_Fraction_', num2str(isoValue), ...
                           '_T', erase(figTime, '.')];

            case 'B'
                figName = ['Far_Field_Instantaneous_Spray_Volume_Fraction_', num2str(isoValue), ...
                           '_T', erase(figTime, '.')];

        end
        
        figSubtitle = [figTime, ' \it{s}'];
        
        fig = plotVolumeField(xLimsData, yLimsData, zLimsData, xInit, yInit, zInit, POD, fieldData, ...
                              fig, figName, geometry, isoValue, cMap, figTitle, figSubtitle, ...
                              xLimsPlot, yLimsPlot, zLimsPlot, figSave);
                           
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
        
        % Format volumeData to Be Human-Readable
        volumeData.positionGrid = [x(:), y(:), z(:)];
        
        for i = 1:nTimes
            volumeData.inst.nParticles{i} = volumeData.inst.nParticles{i}(:);
            volumeData.inst.volFraction{i} = volumeData.inst.volFraction{i}(:);
            volumeData.inst.mass{i} = volumeData.inst.mass{i}(:);
%             volumeData.inst.d32{i} = volumeData.inst.d32{i}(:);
%             volumeData.inst.d30{i} = volumeData.inst.d30{i}(:);
            volumeData.inst.d20{i} = volumeData.inst.d20{i}(:);
            volumeData.inst.d10{i} = volumeData.inst.d10{i}(:);
        end
        clear i;

        volumeData.mean.nParticles = volumeData.mean.nParticles(:);
        volumeData.mean.volFraction = volumeData.mean.volFraction(:);
        volumeData.mean.mass = volumeData.mean.mass(:);
%         volumeData.mean.d32 = volumeData.mean.d32(:);
%         volumeData.mean.d30 = volumeData.mean.d30(:);
        volumeData.mean.d20 = volumeData.mean.d20(:);
        volumeData.mean.d10 = volumeData.mean.d10(:);
        
        volumeData = orderfields(volumeData, {'positionGrid', 'time', 'inst', 'mean'});
        
        % Save Data
        switch format
            
            case 'A'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/nearField'], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/nearField']);
                end
                
            case 'B'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/farField'], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/farField']);
                end
                
        end
        
        switch format
            
            case 'A'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/nearField/', dataID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/nearField/', dataID, '.mat'], ...
                      'caseID', 'dataID', 'volumeData', 'cellSize', 'sampleInterval', 'timePrecision', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/farField/', dataID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/volumeField/', caseID, '/farField/', dataID, '.mat'], ...
                      'caseID', 'dataID', 'volumeData', 'cellSize', 'sampleInterval', 'timePrecision', 'dLims', 'normalise', '-v7.3', '-noCompression');
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