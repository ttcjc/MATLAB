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


%% Lagrangian Particle Tracker v3.1

cloudName = 'kinematicCloud'; % OpenFOAM Cloud Name
ejectionSite = [0.46575, -0.167, 0.006]; % Particle Injection Site Location
normalise = false; % Normalisation of Dimensions

disp('=====================');
disp('Particle Tracker v3.1');
disp('=====================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Updated calls to 'globalPos' to 'positionCartesian'
% v1.2 - Updated to Support Changes to 'timeDirectories.m'
% v2.0 - Rewritten to Support Tracking of Specific Times Instances
% v3.0 - Rewrite, Accommodating New OpenFOAM Data Formats
% v3.1 - 


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
disp('    A: Rear-End Contamination');
disp('    B: Downstream Contaminant Transport');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Region of Interest [A/B]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        formatA = 'A';
        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        formatA = 'B';
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Initialise Lagrangian Data

switch formatA

    case 'A'
        [dataID, LagProps, ~, impactData, volumeData, sampleInterval] = initialiseLagData(saveLocation, caseFolder, caseName, ...
                                                                                          cloudName, false, true, ...
                                                                                          true, timeDirs, deltaT, ...
                                                                                          timePrecision, nProc);

    case 'B'
        [dataID, LagProps, impactData, ~, volumeData, sampleInterval] = initialiseLagData(saveLocation, caseFolder, caseName, ...
                                                                                          cloudName, true, false, ...
                                                                                          true, timeDirs, deltaT, ...
                                                                                          timePrecision, nProc);
        
        % Select Plane of Interest
        planes = fieldnames(impactData);
        
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

        impactData = impactData.(planes{index});
        planePos = erase(planes{index}, '.');
        clear planes;
        
        disp(['Plane of Interest: ', planePos]);
end

disp(' ');
disp(' ');


%% Select Tracking Options

disp('Tracking Options');
disp('-----------------');

disp(' ');

disp('Possible Tracking Methods:');
disp('    A: Particle Ejections at Specified Time Instance');
disp('    B: Particle Impacts at Specified Time Instance');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Tracking Method [A/B]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        formatB = 'A';
        [trackTime, trackTimeIndex] = inputTime('Ejection', impactData.time, impactData.d);
        
        if trackTimeIndex == -1
            continue;
        end

        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        formatB = 'B';
        [trackTime, trackTimeIndex] = inputTime('Impact', impactData.time, impactData.d);

        if trackTimeIndex == -1
            continue;
        end

        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

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

valid = false;
while ~valid
    disp(' ');
    selection = input('Enable Interpolation of Particle Paths? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        interp = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        interp = true;
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Perform Particle Tracking

disp('Particle Tracking');
disp('------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;
evalc('parpool(nProc);');

disp(' ');

disp('    Initialising...');

% Identify Empty Time Instances
i = 1;
while i <= height(volumeData.time)
    
    if isempty(volumeData.d{i})
        
        for j = 1:height(LagProps)
            volumeData.(LagProps{j}){i} = -1;
        end
        clear j;
        
    else
        i = i + 1;
    end
    
end
clear i;

% Remove Unnecessary Data
switch formatB

    % Retain Time Instances > Tracking Time
    case 'A'
        impactData.time = impactData.time((trackTimeIndex + 1):end);
        impactData.timeExact = vertcat(impactData.timeExact{(trackTimeIndex + 1):end});
        
        for i = 1:height(LagProps)
            impactData.(LagProps{i}) = vertcat(impactData.(LagProps{i}){(trackTimeIndex + 1):end});
        end
        clear i;
        
        volumeData.time = volumeData.time((trackTimeIndex - 1):(end - 1));

        for i = 1:height(LagProps)
            volumeData.(LagProps{i})= volumeData.(LagProps{i})((trackTimeIndex - 1):(end - 1));
        end
        clear i;

    % Retain Time Instances < Tracking Time
    case 'B'
        impactData.time = impactData.time(trackTimeIndex);
        impactData.timeExact = impactData.timeExact{trackTimeIndex};
        
        for i = 1:height(LagProps)
            impactData.(LagProps{i}) = impactData.(LagProps{i}){trackTimeIndex};
        end
        clear i;

        volumeData.time = volumeData.time(1:(trackTimeIndex - 1));

        for i = 1:height(LagProps)
            volumeData.(LagProps{i}) = volumeData.(LagProps{i})(1:(trackTimeIndex - 1));
        end
        clear i;

end

% Shift Data Origin
if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
    
    for i = 1:height(volumeData.time)
        
        if volumeData.positionCartesian{i} ~= -1
            volumeData.positionCartesian{i}(:,1) = volumeData.positionCartesian{i}(:,1) + 1.325;
        end

    end
    clear i;
    
    impactData.positionCartesian(:,1) = impactData.positionCartesian(:,1) + 1.325;
end

% Normalise Dimensions
if normalise
    
    if contains(caseName, ["Run_Test", "Windsor"])
        
        for i = 1:height(volumeData.time)
            
            if volumeData.positionCartesian{i} ~= -1
                volumeData.positionCartesian{i} = round((volumeData.positionCartesian{i} / 1.044), ...
                                                        spacePrecision);
            end
            
        end
        clear i;
        
        impactData.positionCartesian = round((impactData.positionCartesian / 1.044), ...
                                             spacePrecision);
    end

end

% Specify Region Boundaries
switch formatA
    
    case 'A'
        
        if contains(caseName, ["Run_Test", "Windsor"])
            xLimsData = [0.31875; 1.26625];
            yLimsData = [-0.3445; 0.3445];
            zLimsData = [0; 0.489];
            
            if normalise
                xLimsData = round((xLimsData / 1.044), spacePrecision);
                yLimsData = round((yLimsData / 1.044), spacePrecision);
                zLimsData = round((zLimsData / 1.044), spacePrecision);
            end

        end
        
    case 'B'
        
        if contains(caseName, ["Run_Test", "Windsor"])
            xLimsData = [0.31875; impactData.positionCartesian(1,1)];
            yLimsData = [-0.5945; 0.5945];
            zLimsData = [0; 0.739];
            
            if normalise
                xLimsData(1) = round((xLimsData(1) / 1.044), spacePrecision);
                yLimsData = round((yLimsData / 1.044), spacePrecision);
                zLimsData = round((zLimsData / 1.044), spacePrecision);
            end
            
        end
        
end

disp(' ');

% Identify Particles of Interest
disp('    Identifying Particles of Interest...');

% Collate Particles of Interest
index = find(((impactData.d * 1e6) >= dLims(1)) & ...
             ((impactData.d * 1e6) <= dLims(2)) & ...
             (impactData.positionCartesian(:,1) >= xLimsData(1)) & ...
             (impactData.positionCartesian(:,1) <= xLimsData(2)) & ...
             (impactData.positionCartesian(:,2) >= yLimsData(1)) & ...
             (impactData.positionCartesian(:,2) <= yLimsData(2)) & ...
             (impactData.positionCartesian(:,3) >= zLimsData(1)) & ...
             (impactData.positionCartesian(:,3) <= zLimsData(2)));

impactData.timeExact = impactData.timeExact(index);
        
for i = 1:height(LagProps)
    impactData.(LagProps{i}) = impactData.(LagProps{i})(index,:);
end
clear i;

% Initialise Progress Bar
wB = waitbar(0, 'Removing Invalid Particles From Volume Data', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, height(volumeData.time));

% Identify Impinging Particles Ejected During Specified Time Instance
switch formatB
    
    case 'A'
        activePreTime = intersect([volumeData.origProcId{1}, volumeData.origId{1}], ...
                                  [impactData.origProcId, impactData.origId], 'rows', 'stable');
        activePostTime = intersect([volumeData.origProcId{2}, volumeData.origId{2}], ...
                                  [impactData.origProcId, impactData.origId], 'rows', 'stable');
                              
        ejections = setdiff(activePostTime, activePreTime, 'rows', 'stable');
        
        index = find(ismember([impactData.origProcId, impactData.origId], ejections, 'rows'));
        
        impactData.timeExact = impactData.timeExact(index);
                
        for i = 1:height(LagProps)
            impactData.(LagProps{i}) = impactData.(LagProps{i})(index,:);
        end
        clear i;
        
        volumeData.time = volumeData.time(2:end);

        for i = 1:height(LagProps)
            volumeData.(LagProps{i})= volumeData.(LagProps{i})(2:end);
        end
        clear i;
        
end

% Remove Invalid Particles From Data
index = cell(height(volumeData.time),1);

particleIDimpact = [impactData.origProcId, impactData.origId];
origProcId = volumeData.origProcId;
origId = volumeData.origId;
parfor i = 1:height(volumeData.time)
    particleIDvolume = [origProcId{i}, origId{i}];

    [~, index{i}] = intersect(particleIDvolume, particleIDimpact, 'rows', 'stable');

    send(dQ, []);
end
clear i particleIDimpact particleIDvolume origProcId origId;

delete(wB);

for i = 1:height(volumeData.time)
    
    for j = 1:height(LagProps)
        volumeData.(LagProps{j}){i} = volumeData.(LagProps{j}){i}(index{i},:);
    end
    clear j;
    
end
clear i;

disp(' ');

% Limit Tracking Count
switch formatA

    case 'A'
        disp(['    ', num2str(height(impactData.positionCartesian)), ' Valid Particles Recorded on Surface'])
    
    case 'B'
        disp(['    ', num2str(height(impactData.positionCartesian)), ' Valid Particles Passed Through Plane of Interest'])

end

intermediateToc = toc;

if height(impactData.timeExact) > 200
    disp('    WARNING: Tracking a Large Number of Particles Will Be Computationally Expensive and Difficult to Visualise');
    disp('             It Is Recommended to Track No More Than 200 Particles');
end

% Specify Particle Count
valid = false;
while ~valid
    disp(' ');
    selection = input('    Limit Number of Particles to Track? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        count = height(impactData.timeExact);

        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        count = inputCount(height(impactData.timeExact));

        if count == -1
            continue;
        end
        
        valid = true;
    else
        disp('        WARNING: Invalid Entry');
    end

end
clear valid;

tic;

% Select a Random Set of Particles to Track
if count ~= height(impactData.timeExact)
    index = sort(randperm(height(impactData.timeExact), count))';
    
    impactData.timeExact = impactData.timeExact(index);

    for i = 1:height(LagProps)
        impactData.(LagProps{i}) = impactData.(LagProps{i})(index,:);
    end
    clear i;
    
end

disp(' ');

% Perform Tracking
disp(['    Tracking ', num2str(count), ' Particles...']);

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Particle Paths', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

trackingData.ID = nan(count,2);
trackingData.d =  nan(count,1);
trackingData.path = cell(count,1);
trackingData.age = trackingData.path;

for i = 1:count
    trackingData.ID(i,:) = [impactData.origProcId(i), impactData.origId(i)];
    trackingData.d(i) = impactData.d(i);
    
    trackingData.path{i} = nan((height(volumeData.time) + 1),3);
    trackingData.age{i} = nan((height(volumeData.time) + 1),1);
    
    for j = 1:height(volumeData.time)
        index = find(ismember([volumeData.origProcId{j}, volumeData.origId{j}], trackingData.ID(i,:), 'rows'));
        
        if ~isempty(index)
            trackingData.path{i}((j + 1),:) = volumeData.positionCartesian{j}(index,:);
            
            if isnan(trackingData.age{i}(j))
                trackingData.age{i}(j + 1) = deltaT;
            else
                trackingData.age{i}(j + 1) = trackingData.age{i}(j) + deltaT;
            end
            
        end
        
    end
    clear j;
    
    % Shift Initial Position to Ejection Site
    index = find(isnan(trackingData.age{i}) == false, 1, 'first');
    
    trackingData.path{i}((index - 1),:) = ejectionSite;
    trackingData.age{i}(index - 1) = 0;
    
    % Shift Final Position to Impact Site
    index = find(isnan(trackingData.age{i}) == false, 1, 'last');

    switch formatA

        case 'A'
            trackingData.path{i}((index + 1),1) = xDims(2) + 1e-6;
            trackingData.path{i}((index + 1),(2:3)) = impactData.positionCartesian(i,(2:3));

        case 'B'
            trackingData.path{i}((index + 1),:) = impactData.positionCartesian(i,:);

    end

    trackingData.age{i}(index + 1) = trackingData.age{i}(index) + deltaT;
    
    waitbar((i / count), wB);
end
clear i;

delete(wB);

% Interpolate Particle Trajectories
if interp
    % Initialise Progress Bar
    wB = waitbar(0, 'Interpolating Particle Trajectories', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    
    for i = 1:count
        index = find(isnan(trackingData.age{i}) == false);
        
        trackingData.path{i} = trackingData.path{i}(index,:);
        trackingData.age{i} = trackingData.age{i}(index);
    
        trackingData.path{i} = interparc((2 * height(trackingData.path{i})), ...
                                         trackingData.path{i}(:,1), ...
                                         trackingData.path{i}(:,2), ...
                                         trackingData.path{i}(:,3), 'spline');
        
        waitbar((i / count), wB);
    end
    clear i;
    
end

delete(wB);

evalc('delete(gcp(''nocreate''));');
executionTime = intermediateToc + toc;

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
    selection = input('Plot Particle Paths? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotPaths = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotPaths = true;
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

if plotPaths
    colourVar = {'Diameter', 'Age'};
    
    valid = false;
    while ~valid
        disp(' ');
        [index, valid] = listdlg('listSize', [300, 300], ...
                                 'selectionMode', 'single', ...
                                 'name', 'Select Variable Used for Particle Colouration', ...
                                 'listString', colourVar);

        if ~valid
            disp('WARNING: No Mapping Variable Selected');
        end

    end
    clear valid;

    colourVar = colourVar{index};
    
    if strcmp(colourVar, 'Diameter')
        minColourVar = dLims(1);
        maxColourVar = dLims(2);
    elseif strcmp(colourVar, 'Age')
        minColourVar = deltaT;
        maxColourVar = max(cellfun(@max, trackingData.age));
    end

    disp(['Variable of Interest: ', colourVar]);
end

disp(' ');
disp(' ');


%% Present Particle Paths

disp('Particle Path Presentation');
disp('---------------------------');

disp(' ');

if plotPaths

    switch formatA

        case 'A'

            switch formatB

                case 'A'
                    figName = ['Rear_End_Soiling_Time_Eject_T', erase(num2str(trackTime), '.'), ...
                               '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];

                case 'B'
                    figName = ['Rear_End_Soiling_Time_Impact_T', erase(num2str(trackTime), '.'), ...
                               '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];

            end

        case 'B'

            switch formatB

                case 'A'
                    figName = ['Downstream_Spray_Time_Eject_T', erase(num2str(trackTime), '.'), ...
                               '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];

                case 'B'
                    figName = ['Downstream_Spray_Time_Impact_T', erase(num2str(trackTime), '.'), ...
                               '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];

            end

    end

    cMap = viridis(32);
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    figSubtitle = ' ';
    xLimsPlot = xLimsData;
    yLimsPlot = yLimsData;
    zLimsPlot = zLimsData;
    
    fig = plotParticlePaths(trackingData, fig, figName, geometry, cMap, colourVar, ...
                            minColourVar, maxColourVar, figTitle, figSubtitle, ...
                            xLimsPlot, yLimsPlot, zLimsPlot);
else
    disp('    Skipping Volume Field Presentation');
end


%% Local Functions

function [time, index] = inputTime(type, timeList, dList)

    time = str2double(input(['    Input Desired ', type, ' Time [s]: '], 's'));

    index = find(ismember(timeList, time));
    
    if isempty(time)
        disp('        WARNING: Invalid Entry');
        index = -1;
    elseif isempty(dList{index})
        disp(['        WARNING: No ', type, 's Recorded During Specified Time Instance']);
        index = -1;
    end
    
end


function D = inputD(type)

    D = str2double(input(['    ', type, ' Diameter of Interest [', char(956), 'm]: '], 's'));
    
    if isnan(D) || length(D) > 1 || D < 1
        disp('        WARNING: Invalid Entry');
        D = -1;
    end
    
end


function count = inputCount(nParticles)
    
    count = str2double(input(['        Input Number of Particles to Track [1-', num2str(nParticles), ']: '], 's'));
    
    if isnan(count) || count <= 0 || count >= nParticles
        disp('            WARNING: Invalid Entry');
        count = -1;
    end
    
end