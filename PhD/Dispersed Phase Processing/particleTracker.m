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


%% Lagrangian Particle Tracker v4.0

cloudName = 'kinematicCloud'; % OpenFOAM Cloud Name

figSave = false; % Save .fig File(s)

normalise = false; % Normalisation of Dimensions

disp('=====================');
disp('Particle Tracker v4.0');
disp('=====================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Updated calls to 'globalPos' to 'positionCartesian'
% v1.2 - Updated to Support Changes to 'timeDirectories.m'
% v2.0 - Rewritten to Support Tracking of Specific Times Instances
% v3.0 - Rewrite, Accommodating New OpenFOAM Data Formats
% v3.1 - Added Support for Full-Scale Windsor Model Simulations
% v4.0 - Enabled Tracking of Particles That Impact User-Specified Region


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
disp('    A: Rear-End Contamination');
disp('    B: Downstream Contaminant Transport');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Region of Interest [A/B]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        format1 = '1A';
        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        format1 = '1B';
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

switch format1

    case '1A'
        [dataID, LagProps, ~, impactData, volumeData, sampleInterval] = initialiseLagData(saveLocation, caseFolder, caseID, ...
                                                                                          cloudName, false, true, ...
                                                                                          true, timeDirs, deltaT, ...
                                                                                          timePrecision, maxNumCompThreads);

    case '1B'
        [dataID, LagProps, impactData, ~, volumeData, sampleInterval] = initialiseLagData(saveLocation, caseFolder, caseID, ...
                                                                                          cloudName, true, false, ...
                                                                                          true, timeDirs, deltaT, ...
                                                                                          timePrecision, maxNumCompThreads);
        
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
disp('    A: Ejections at Specified Time Instance');
disp('    B: Impacts at Specified Time Instance');
disp('    C: Impacts in Specified Region');

figHold = figHold + fig;

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Tracking Method [A/B/C]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        format2 = '2A';
        
        [trackTime, trackTimeIndex] = inputTime('Ejection', impactData.time, impactData.d);
        
        if trackTimeIndex == -1
            continue;
        end

        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        format2 = '2B';
        
        [trackTime, trackTimeIndex] = inputTime('Impact', impactData.time, impactData.d);

        if trackTimeIndex == -1
            continue;
        end

        valid = true;
    elseif selection == 'c' | selection == 'C' %#ok<OR2>
        format2 = '2C';
        
        fig = figHold + 1;
        impactPosition = inputPosition(fig, unique(vertcat(impactData.positionCartesian{:}), 'rows'));

        if impactPosition == -1
            continue;
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

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

evalc('parpool(''threads'');'); % evalc('parpool(nProc);');

%%%%

disp(' ');

disp('    Initialising...');

% Identify Empty Time Instances Not Related to Case Initialisation
emptyTimes = cellfun(@isempty, volumeData.origId);
firstValidTime = find(emptyTimes == false, 1, 'first');

suspectTimes = find(emptyTimes(firstValidTime:end) == true) + (firstValidTime - 1);

if ~isempty(suspectTimes)
    
    for i = 1:height(suspectTimes)
        disp(['        WARNING: Time ''', num2str(volumeData.time(suspectTimes(i))), ''' Is Unexpectedly Empty']);
    end
    clear i;
    
end

% Remove Unnecessary Data
impactData = rmfield(impactData, 'time');

switch format2

    % Retain Time Instances > Tracking Time
    case '2A'
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
    case '2B'
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

    % Retain All Time Instances
    case '2C'
        impactData.timeExact = vertcat(impactData.timeExact{:});
        
        for i = 1:height(LagProps)
            impactData.(LagProps{i}) = vertcat(impactData.(LagProps{i}){:});
        end
        clear i;
        
end

nTimes = height(volumeData.time);

% Shift Data Origin
if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
    
    for i = 1:nTimes
        
        if ~isempty(volumeData.positionCartesian{i})
            volumeData.positionCartesian{i}(:,1) = volumeData.positionCartesian{i}(:,1) + 1.325;
        end
        
    end
    clear i;
    
    impactData.positionCartesian(:,1) = impactData.positionCartesian(:,1) + 1.325;
    
    switch format2
        
        case '2C'
            impactPosition(1) = impactPosition(1) + 1.325;
            
    end
    
end

% Normalise Coordinate System
if normalise
    disp('        Normalising Coordinate System');
    
    % Initialise Progress Bar
    wB = waitbar(0, 'Normalising Coordinate System', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';

    % Perform Normalisation
    for i = 1:nTimes

        if ~isempty(volumeData.positionCartesian{i})
            volumeData.positionCartesian{i}  = round((volumeData.positionCartesian{i} / normLength), ...
                                                  spacePrecision);
        end

        % Update Waitbar
        waitbar((i / nTimes), wB);
    end
    clear i;
    
    impactData.positionCartesian = round((impactData.positionCartesian / normLength), ...
                                         spacePrecision);

    delete(wB);
end

% Specify Region Boundaries
switch format2
    
    case '2C'
        
        if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
            regionSize = 16e-3;
        else
            regionSize = 64e-3;
        end
        
        xLimsData = impactPosition(1);
        yLimsData = [(impactPosition(2) - regionSize); (impactPosition(2) + regionSize)];
        zLimsData = [(impactPosition(3) - regionSize); (impactPosition(3) + regionSize)];
        
    otherwise
        
        switch format1

            case '1A'

                if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
                    xLimsData = impactData.positionCartesian(1,1);
                    yLimsData = [-0.4176; 0.4176];
                    zLimsData = [0; 0.4176];
                else
                    xLimsData = impactData.positionCartesian(1,1);
                    yLimsData = [-1.6704; 1.6704];
                    zLimsData = [0; 1.6704];
                end

            case '1B'

                if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
                    xLimsData = impactData.positionCartesian(1,1);
                    yLimsData = [-0.6264; 0.6264];
                    zLimsData = [0; 0.6264];
                else
                    xLimsData = impactData.positionCartesian(1,1);
                    yLimsData = [-2.088; 2.088];
                    zLimsData = [0; 2.088];
                end

        end
        
end
        
if normalise
    xLimsData = round((xLimsData / normLength), spacePrecision);
    yLimsData = round((yLimsData / normLength), spacePrecision);
    zLimsData = round((zLimsData / normLength), spacePrecision);
end

disp(' ');

% Collate Particles of Interest
disp('    Collating Particles of Interest...');

% Collate Particles of Interest
index = find(((impactData.d * 1e6) >= dLims(1)) & ...
             ((impactData.d * 1e6) <= dLims(2)) & ...
             (impactData.positionCartesian(:,2) >= yLimsData(1)) & ...
             (impactData.positionCartesian(:,2) <= yLimsData(2)) & ...
             (impactData.positionCartesian(:,3) >= zLimsData(1)) & ...
             (impactData.positionCartesian(:,3) <= zLimsData(2)));

impactData.timeExact = impactData.timeExact(index);
        
for i = 1:height(LagProps)
    impactData.(LagProps{i}) = impactData.(LagProps{i})(index,:);
end
clear i;

% Temporal Filtering
switch format2
    
    % Identify Impinging Particles Ejected During Specified Time Instance
    case '2A'
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
        
    % Remove Particles That Existed During First Time Instance
    otherwise
        activeFirstTime = intersect([volumeData.origProcId{1}, volumeData.origId{1}], ...
                                    [impactData.origProcId, impactData.origId], 'rows', 'stable');
        
        index = find(ismember([impactData.origProcId, impactData.origId], activeFirstTime, 'rows'));
        
        impactData.timeExact(index) = [];
                
        for i = 1:height(LagProps)
            impactData.(LagProps{i})(index,:) = [];
        end
        clear i;

end

volumeData.time = volumeData.time(2:end);

nTimes = height(volumeData.time);

for i = 1:height(LagProps)
    volumeData.(LagProps{i})= volumeData.(LagProps{i})(2:end);
end
clear i;

disp(' ');

% Limit Tracking Count
switch format1

    case '1A'
        disp(['    ', num2str(height(impactData.positionCartesian)), ' Valid Particles Recorded on Surface'])
    
    case '1B'
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

% Remove Invalid Particles From Data
disp('    Removing Invalid Particles From Data...');

% Initialise Progress Bar
wB = waitbar(0, 'Removing Invalid Particles From Data', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, nTimes);

% Perform Removal
index = cell(nTimes,1);

particleIDimpact = [impactData.origProcId, impactData.origId];
origProcId = volumeData.origProcId;
origId = volumeData.origId;
parfor i = 1:nTimes
    particleIDvolume = [origProcId{i}, origId{i}];

    [~, index{i}] = intersect(particleIDvolume, particleIDimpact, 'rows', 'stable');

    % Remove Unnecessary Data
    origProcId{i} = [];
    origId{i} = [];
    
    % Update Waitbar
    send(dQ, []);
end
clear i particleIDimpact particleIDvolume origProcId origId;

delete(wB);

for i = 1:nTimes
    
    for j = 1:height(LagProps)
        volumeData.(LagProps{j}){i} = volumeData.(LagProps{j}){i}(index{i},:);
    end
    clear j;
    
end
clear i;

disp(' ');

% Perform Tracking
disp(['    Tracking ', num2str(count), ' Particles...']);

if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
    ejectionSite = [0.46575, -0.167, 0.006];
else
    ejectionSite = [];
end

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Particle Paths', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Perform Tracking
trackingData.ID = nan(count,2);
trackingData.d =  nan(count,1);
trackingData.path = cell(count,1);
trackingData.age = trackingData.path;

for i = 1:count
    trackingData.ID(i,:) = [impactData.origProcId(i), impactData.origId(i)];
    trackingData.d(i) = impactData.d(i) * 1e6;
    
    trackingData.path{i} = nan((nTimes + 1),3);
    trackingData.age{i} = nan((nTimes + 1),1);
    
    for j = 1:nTimes
        index = find(ismember([volumeData.origProcId{j}, volumeData.origId{j}], trackingData.ID(i,:), 'rows'));
        
        if ~isempty(index)
            
            switch format1
                
                case '1B'
                    
                    if volumeData.positionCartesian{j}(index,1) > xLimsData
                        continue;
                    end
                    
            end
            
            trackingData.path{i}((j + 1),:) = volumeData.positionCartesian{j}(index,:);
            
            if isnan(trackingData.age{i}(j))
                trackingData.age{i}(j) = 0;
            end
            
            trackingData.age{i}(j + 1) = trackingData.age{i}(j) + deltaT;
        end
        
    end
    clear j;
    
    % Shift Initial Position to Ejection Site
    if ~isempty(ejectionSite)
        trackingData.path{i}((trackingData.age{i} == 0),:) = ejectionSite;
    end
    
    % Shift Final Position to Impact Site
    index = find(isnan(trackingData.age{i}) == false, 1, 'last');
    
    switch format1

        case '1A'
            trackingData.path{i}((index + 1),1) = xDims(2) + 1e-3; % Offset Particles From Base for Better Visibility
            trackingData.path{i}((index + 1),(2:3)) = impactData.positionCartesian(i,(2:3));

        case '1B'
            trackingData.path{i}((index + 1),:) = impactData.positionCartesian(i,:);

    end

    % Calculate Final Particle Age
    trackingData.age{i}(index + 1) = trackingData.age{i}(index) + deltaT;
    
    % Remove Unnecessary Time Instances 
    index = find(isnan(trackingData.path{i}(:,1)) == false);

    trackingData.path{i} = trackingData.path{i}(index,:);
    trackingData.age{i} = trackingData.age{i}(index);
    
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
        trackingData.path{i} = interparc((4 * height(trackingData.path{i})), ...
                                         trackingData.path{i}(:,1), ...
                                         trackingData.path{i}(:,2), ...
                                         trackingData.path{i}(:,3), 'spline');
                                     
        trackingData.age{i} = interparc((4 * height(trackingData.age{1})), ...
                                        (1:height(trackingData.age{1})), ...
                                        trackingData.age{1}, 'spline');
        trackingData.age{i} = trackingData.age{i}(:,2);
        
        
        waitbar((i / count), wB);
    end
    clear i;
    
end

delete(wB);

%%%%

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
        minColourVar = 0;
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

    switch format1

        case '1A'

            switch format2

                case '2A'
                    figName = ['Rear_End_Soiling_Time_Eject_T', erase(num2str(trackTime), '.'), ...
                               '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];
                    
                case '2B'
                    figName = ['Rear_End_Soiling_Time_Impact_T', erase(num2str(trackTime), '.'), ...
                               '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];
                    
                case '2C'
                    posX = round(impactPosition(1), 1);
                    posY = round(impactPosition(2), 1);
                    posZ = round(impactPosition(3), 1);
                    
                    if posX >= 0
                        posX = ['P', erase(num2str(posX), '.')];
                    else
                        posX = ['N', erase(num2str(posX), [".", "-"])];
                    end
                    
                    if posY >= 0
                        posY = ['P', erase(num2str(posY), '.')];
                    else
                        posY = ['N', erase(num2str(posY), [".", "-"])];
                    end
                    
                    if posZ >= 0
                        posZ = ['P', erase(num2str(posZ), '.')];
                    else
                        posZ = ['N', erase(num2str(posZ), [".", "-"])];
                    end
                    
                    figName = ['Rear_End_Soiling_Impacts_x', posX,'_y', posY, '_z', posZ];

            end

        case '1B'

            switch format2

                case '2A'
                    figName = ['Downstream_Spray_Time_Eject_T', erase(num2str(trackTime), '.'), ...
                               '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];
                    
                case '2B'
                    figName = ['Downstream_Spray_Time_Impact_T', erase(num2str(trackTime), '.'), ...
                               '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];
                    
                case '2C'
                    posX = round(impactPosition(1), 1);
                    posY = round(impactPosition(2), 1);
                    posZ = round(impactPosition(3), 1);
                    
                    if posX >= 0
                        posX = ['P', erase(num2str(posX), '.')];
                    else
                        posX = ['N', erase(num2str(posX), [".", "-"])];
                    end
                    
                    if posY >= 0
                        posY = ['P', erase(num2str(posY), '.')];
                    else
                        posY = ['N', erase(num2str(posY), [".", "-"])];
                    end
                    
                    if posZ >= 0
                        posZ = ['P', erase(num2str(posZ), '.')];
                    else
                        posZ = ['N', erase(num2str(posZ), [".", "-"])];
                    end
                    
                    figName = ['Downstream_Spray_Impacts_x', posX,'_y', posY, '_z', posZ];

            end

    end

    cMap = viridis(32);
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    figSubtitle = ' ';
    
    switch format1

        case '1A'

            if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
                xLimsPlot = [0.31875; 1.43075];
                zLimsPlot = [0; 0.4176];
            else
                xLimsPlot = [1.275; 5.723];
                yLimsPlot = [-1.6704; 1.6704];
                zLimsPlot = [0; 1.6704];
            end
            
            if normalise
                xLimsPlot = round((xLimsPlot / normLength), spacePrecision);
                yLimsPlot = round((yLimsPlot / normLength), spacePrecision);
                zLimsPlot = round((zLimsPlot / normLength), spacePrecision);
            end

        case '1B'

            if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
                xLimsPlot = [0.31875; impactData.positionCartesian(1,1)];
                yLimsPlot = [-0.6264; 0.6264];
                zLimsPlot = [0; 0.6264];
            else
                xLimsPlot = [1.275; impactData.positionCartesian(1,1)];
                yLimsPlot = [-2.088; 2.088];
                zLimsPlot = [0; 2.088];
            end
            
            if normalise
                xLimsPlot(1) = round((xLimsPlot(1) / normLength), spacePrecision);
                yLimsPlot = round((yLimsPlot / normLength), spacePrecision);
                zLimsPlot = round((zLimsPlot / normLength), spacePrecision);
            end

    end
    
    fig = plotParticlePaths(trackingData, fig, figName, geometry, cMap, colourVar, ...
                            minColourVar, maxColourVar, figTitle, figSubtitle, ...
                            xLimsPlot, yLimsPlot, zLimsPlot, figSave);
else
    disp('    Skipping Volume Field Presentation');
end


%% Local Functions

function [time, index] = inputTime(type, timeList, dList)

    time = single(str2double(input(['    Input Desired ', type, ' Time [s]: '], 's')));

    index = find(ismember(timeList, time));
    
    if isempty(index)
        disp('        WARNING: Invalid Entry');
        
        time = -1; index = -1;
    elseif index == 1
        disp('        WARNING: Particle Ejections During First Time Instance Cannot Be Verified');
        
        time = -1; index = -1;
    elseif isempty(dList{index})
        disp(['        WARNING: No ', type, 's Recorded During Specified Time Instance']);
        
        time = -1; index = -1;
    end
    
end


function [pos] = inputPosition(fig, impactList)
    
    yBoundary = ceil(max(abs(impactList(:,2))) * 2) / 2;
    zBoundary = ceil(max(abs(impactList(:,3))) * 2) / 2;
    
    disp('    Select an Impact Point [Right-Click]:')
    
    set(figure(fig), 'name', 'Particle Impacts', 'color', [1, 1, 1], 'units', 'pixels');
    set(gca, 'dataAspectRatio', [1, 1, 1], 'lineWidth', 4, 'fontName', 'LM Mono 12', ...
             'fontSize', 18, 'layer', 'top');
    hold on;
    
    cMap = viridis(3); cMap = cMap(1,:);
    
    scatter(impactList(:,2), impactList(:,3), 5, cMap, 'filled');
    
    xlim([-yBoundary, yBoundary]);
    ylim([0, zBoundary]);
    tickData = -yBoundary:((yBoundary - -yBoundary) / 4):yBoundary;
    xticks(tickData(2:(end-1)));
    tickData = 0:(zBoundary / 4):zBoundary;
    yticks(tickData(2:(end-1)));
    hold off;
    
    try
        [yPos, zPos] = getpts;
        
        if numel(yPos) > 1
            disp('        Warning: Multiple Impact Positions Selected');
            pos = -1;
            close(fig);
        elseif yPos <= min(impactList(:,2)) || yPos >= max (impactList(:,2)) || ...
               zPos <= min(impactList(:,3)) || zPos >= max (impactList(:,3))
            disp('        Warning: Selected Position Must Be Within Impact Area');
            
            pos = -1;
            close(fig);
        else
            pos = [impactList(1,1), yPos, zPos];
            close(fig);
        end
        
    catch
        disp('        Warning: A Valid Impact Position Must Be Selected');
        
        pos = -1;
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
