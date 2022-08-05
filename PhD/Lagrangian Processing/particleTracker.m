%% Lagrangian Particle Tracker v3.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

saveLocation = '/mnt/Processing/Data';
% saveLocation = '~/Data';

normalise = true; % Normalisation of Dimensions

cloudName = 'kinematicCloud'; % OpenFOAM Cloud Name

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp('=====================');
disp('Particle Tracker v3.0');
disp('=====================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Updated calls to 'globalPos' to 'positionCartesian'
% v1.2 - Updated to Support Changes to 'timeDirectories.m'
% v2.0 - Rewritten to Support Tracking of Specific Times Instances
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

% if normalise
%     dataID = [dataID, '_Norm'];
% end

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
        formatB = 'C';
        trackingTime = inputTime('Ejection', impactData.time, impactData.d);
        
        if trackingTime == -1
            continue;
        end

        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        formatB = 'D';
        trackingTime = inputTime('Impact', impactData.time, impactData.d);

        if trackingTime == -1
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

% if normalise
%     dataID = insertBefore(dataID, '_Norm', ['_D', num2str(dLims(1)), '_D', num2str(dLims(2))]);
% else
%     dataID = [dataID, '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];
% end

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
        
    else
        i = i + 1;
    end
    
end
clear i;

% Remove Unnecessary Data
switch formatB

    % Retain Time Instances > Tracking Time
    case 'C'
        impactData.time = impactData.time((trackingTime + 1):end);
        impactData.timeExact = vertcat(impactData.timeExact{(trackingTime + 1):end});
        
        for i = 1:height(LagProps)
            impactData.(LagProps{i}) = vertcat(impactData.(LagProps{i}){(trackingTime + 1):end});
        end
        
        volumeData.time = volumeData.time((trackingTime - 1):(end - 1));

        for i = 1:height(LagProps)
            volumeData.(LagProps{i})= volumeData.(LagProps{i})((trackingTime - 1):(end - 1));
        end

    % Retain Time Instances < Tracking Time
    case 'D'
        impactData.time = impactData.time(trackingTime);
        impactData.timeExact = impactData.timeExact{trackingTime};
        
        for i = 1:height(LagProps)
            impactData.(LagProps{i}) = impactData.(LagProps{i}){trackingTime};
        end

        volumeData.time = volumeData.time(1:(trackingTime - 1));

        for i = 1:height(LagProps)
            volumeData.(LagProps{i}) = volumeData.(LagProps{i})(1:(trackingTime - 1));
        end

end

% Shift Data Origin
if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
    
    for i = 1:height(volumeData.time)
        
        if volumeData.positionCartesian{i} ~= -1
            volumeData.positionCartesian{i}(:,1) = volumeData.positionCartesian{i}(:,1) + 1.325;
        end

    end
    
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
        
        impactData.positionCartesian = round((impactData.positionCartesian / 1.044), ...
                                             spacePrecision);
    end

end

% Specify Region Boundaries
switch formatA
    
    case 'A'
        
        if contains(caseName, ["Run_Test", "Windsor"])
            xLimsData = [0.31875; 1.38325];
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
            zLimsData = [0; 0.489];
            
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

% Initialise Progress Bar
wB = waitbar(0, 'Removing Invalid Particles From Volume Data', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, height(volumeData.time));

% Identify Impinging Particles Ejected During Specified Time Instance
switch formatB
    
    case 'C'
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
        
        volumeData.time = volumeData.time(2:end);

        for i = 1:height(LagProps)
            volumeData.(LagProps{i})= volumeData.(LagProps{i})(2:end);
        end
        
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
clear particleIDimpact particleIDvolume origProcId origId;

delete(wB);

for i = 1:height(volumeData.time)
    
    for j = 1:height(LagProps)
        volumeData.(LagProps{j}){i} = volumeData.(LagProps{j}){i}(index{i},:);
    end
    
end

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

% Specify Sampling Frequency
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
    
end

% Perform Tracking
disp(['    Tracking ', num2str(count), ' Particles...']);

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Particle Paths', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

trackingData.ID = nan(count,2);
trackingData.d =  nan(count,1);
trackingData.path = cell(count,height(volumeData.time));
trackingData.age = nan(count,height(volumeData.time));

for i = 1:count
    trackingData.ID(i,:) = [impactData.origProcId(i), impactData.origId(i)];
    trackingData.d(i) = impactData.d(i);
    
    for j = 1:height(volumeData.time)
        index = find(ismember([volumeData.origProcId{j}, volumeData.origId{j}], trackingData.ID(i,:), 'rows'));
        
        if ~isempty(index)
            trackingData.path{i,j} = volumeData.positionCartesian{j}(index,:);
            
            if j == 1 || isnan(trackingData.age(i,(j - 1)))
                trackingData.age(i,j) = deltaT;
            else
                trackingData.age(i,j) = trackingData.age(i,(j - 1)) + deltaT;
            end
            
        end
        
    end
    
    waitbar((i / count), wB);
end

% Add Initial Injection Location
% Add Final Impact Location

delete(wB);

% Make Pretty Pictures

disp(' ');
disp(' ');


%% Local Functions

function time = inputTime(type, timeList, dList)

    time = str2double(input(['    Input Desired ', type, ' Time [s]: '], 's'));

    time = find(ismember(timeList, time));
    
    if isempty(time)
        disp('        WARNING: Invalid Entry');
        time = -1;
    elseif isempty(dList{time})
        disp(['        WARNING: No ', type, 's Recorded During Specified Time Instance']);
        time = -1;
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