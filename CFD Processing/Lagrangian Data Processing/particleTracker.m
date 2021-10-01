%% Lagrangian Particle Tracker v2.0

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('=====================');
disp ('Particle Tracker v2.0');
disp ('=====================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Updated calls to 'globalPos' to 'positionCartesian'
% v1.2 - Updated to Support Changes to 'timeDirectories.m'
% v2.0 - Rewritten to Support Tracking of Specific Times Instances


%% Case Initialisation

[caseFolder, xDims, yDims, zDims, timeDirs, deltaT, geometry] = initialiseCase('global');

disp(' ');
disp(' ');


%% Lagrangian Data Acquisition

disp('LAGRANGIAN DATA ACQUISITION');
disp('---------------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Load Saved Particle Data? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        disp(' ');
        [particleData, particleProps] = lagrangianData(caseFolder, timeDirs);
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        [fileName, filePath] = uigetfile('/mnt/Processing/Data/Numerical/MATLAB/particleData/*.*', ...
                                         'Select Lagrangian Data');

        if contains(filePath, '/MATLAB/particleData/') 
            variableInfo = who('-file', [filePath, fileName]);
            
            if ismember('particleDataGlobal', variableInfo)
                disp(['    Loading: ', fileName]);
                load([filePath, fileName]);
                particleData = particleDataGlobal;
                particleProps = particlePropsGlobal;
                clearvars particleDataGlobal particlePropsGlobal
                valid = true;
            else
                disp('    WARNING: Global Lagrangian Data Not Available in Specified File');
                clearvars fileName filePath;
            end
            
        else
            disp('    WARNING: Invalid File Selection');
            clearvars fileName filePath;
        end

    else
        disp('    WARNING: Invalid Entry');
    end

end

clearvars variableInfo;

disp(' ');
disp(' ');


%% Tracking Options

disp('TRACKING OPTIONS');
disp('----------------');

disp(' ');

disp('Possible Tracking Methods:');
disp('    A: Base Contamination');
disp('    B: Full-Body Contamination (NYI)');
disp('    C: Far-Field Particle Transport');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Mapping Location [A/B/C]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        method = 'A';
        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        method = 'B';
        valid = true;
    elseif selection == 'c' | selection == 'C' %#ok<OR2>
        method = 'C';
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end

switch method
    
    case {'A', 'B'}
        disp(' ');
        
        disp('Lagrangian Data Available for the Following Time Instances:');
        
        for i = 1:height(particleData.time)
            disp(['        T = ', num2str(str2double(particleData.time{i,1}), '%.4f'), ' s']);
        end
        
        valid = false;
        while ~valid
            disp(' ');
            selection = input('Track Impacts at Specific Time Instance? [y/n]: ', 's');

            if selection == 'n' | selection == 'N' %#ok<OR2>
                valid = true;
            elseif selection == 'y' | selection == 'Y' %#ok<OR2>
                impactTime = inputTime(particleData.time);
                
                if isempty(impactTime)
                    disp('        WARNING: No Lagrangian Data Available for Selected Time Instance');
                elseif impactTime == 1
                    disp('        WARNING: Unable to Track Particles Impacts for First Time Instance');
                else
                    valid = true;
                end
                
            else
                disp('    WARNING: Invalid Entry');
            end

        end
        
end

valid = false;
while ~valid
    disp(' ');
    selection = input('Filter Particle Diameters? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        minD = floor(min(particleData.d{end,1}) * 1e6);
        maxD = ceil(max(particleData.d{end,1}) * 1e6);
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        minD = inputD('Minimum');
        maxD = inputD('Maximum');

        if maxD < floor(min(particleData.d{end,1}) * 1e6) || minD > ceil(max(particleData.d{end,1}) * 1e6)
            disp('        WARNING: No Lagrangian Data in Selected Data Range');
        elseif maxD < minD
            disp('        WARNING: Invalid Entry');
        else
            valid = true;
        end

    else
        disp('    WARNING: Invalid Entry');
    end

end

valid = false;
while ~valid
    disp(' ');
    selection = input('Enable Position Interpolation? [y/n]: ', 's');

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


disp(' ');
disp(' ');

% Temporary Implementation Control
if method == 'B' 
    error('Mapping Location Not Yet Implemented');
end


%% Particle Tracking

disp('PARTICLE TRACKING');
disp('-----------------');

disp(' ');

disp('Initialising... ');

% Shift Data Origin
if contains(caseFolder, 'Upstream')
    xDims = xDims + 1.325;
    
    for i = 1:height(particleData.time)
        particleData.positionCartesian{i,1}(:,1) = particleData.positionCartesian{i,1}(:,1) + 1.325;
    end
    
else
    % Add Support for Future Configurations
end

% Set and Normalise Dimensions
xPre = max(width(extractAfter(num2str(xDims(1), 8), '.')), width(extractAfter(num2str(xDims(2), 8), '.')));
yPre = max(width(extractAfter(num2str(yDims(1), 8), '.')), width(extractAfter(num2str(yDims(2), 8), '.')));
zPre = max(width(extractAfter(num2str(zDims(1), 8), '.')), width(extractAfter(num2str(zDims(2), 8), '.')));

if contains(caseFolder, 'Test_Block') || contains(caseFolder, 'Windsor')    
    xDims = round(xDims / 1.044, xPre);
    yDims = round(yDims / 1.044, yPre);
    zDims = round(zDims / 1.044, zPre);
    
    switch method
        
        case 'A'
            xLims = round([0.31875; 1.075] / 1.044, xPre);
            yLims = round([-0.3945; 0.3945] / 1.044, yPre);
            zLims = round([0; 0.539] / 1.044, zPre);
            
        case 'B'
            % Add Support for Full-Body Tracking
            
        case 'C'
            xLims = round([0.31875; 2.573] / 1.044, xPre);
            yLims = round([-0.3945; 0.3945] / 1.044, yPre);
            zLims = round([0; 0.539] / 1.044, zPre);
            
    end
    
    for i = 1:height(particleData.time)
        particleData.positionCartesian{i,1}(:,1) = round(particleData.positionCartesian{i,1}(:,1) / 1.044, xPre);
        particleData.positionCartesian{i,1}(:,2) = round(particleData.positionCartesian{i,1}(:,2) / 1.044, yPre);
        particleData.positionCartesian{i,1}(:,3) = round(particleData.positionCartesian{i,1}(:,3) / 1.044, zPre);
    end
    
else
    % Add Support for Future Geometries
end

% Identify Model Boundaries
switch method
    
    case 'A'
        part = fieldnames(geometry);
        for i = 1:height(part)

            if round(max(geometry.(part{i,1}).vertices(:,1)), xPre) == xDims(2)
                break
            end

            if i == height(part)
                disp(' ');
                error('Mismatch Between ''xDims'' and Geometry Bounding Box')
            end

        end

        geoPoints = unique(geometry.(part{i,1}).vertices, 'stable', 'rows');
        index = find(round(geoPoints(:,1), xPre) == xDims(2));

        xDimsBase = xDims(2);
        yDimsBase = round([min(geoPoints(index,2)); max(geoPoints(index,2))], yPre);
        zDimsBase = round([min(geoPoints(index,3)); max(geoPoints(index,3))], zPre);
        
        
    case 'B'
        % Add Support for Full-Body Mapping

end

disp(' ');

% Identify All Impinging Particles
disp('Identifying Particles of Interest...');

switch method
    
    case 'A'
        index = find(~particleData.active{end,1} & ...
                     (particleData.d{end,1} * 1e6) >= minD & ...
                     (particleData.d{end,1} * 1e6) <= maxD & ...
                     particleData.positionCartesian{end,1}(:,1) == xDimsBase & ...
                     particleData.positionCartesian{end,1}(:,2) >= yDimsBase(1) & ...
                     particleData.positionCartesian{end,1}(:,2) <= yDimsBase(2) & ...
                     particleData.positionCartesian{end,1}(:,3) >= zDimsBase(1) & ...
                     particleData.positionCartesian{end,1}(:,3) <= zDimsBase(2));
                 
    case 'B'
        % Add Support for Full-Body Mapping
            
    case 'C'
        index = find(~particleData.active{end,1} & ...
                     (particleData.d{end,1} * 1e6) >= minD & ...
                     (particleData.d{end,1} * 1e6) <= maxD & ...
                     particleData.positionCartesian{end,1}(:,1) >= xDims(2) +  1);
                 
end

% Collate All Impinging Particles
for i = 1:height(particleProps)
    prop = particleProps{i,1};
    contaminantData.(prop) = particleData.(prop){end,1}(index,:);
end

% Identify Particles of Interest
if exist('impactTime', 'var')
    index = find(~particleData.active{impactTime,1});
    impactsCurrent = intersect(horzcat(particleData.origId{impactTime,1}(index,1), particleData.origProcId{impactTime,1}(index,1)), ...
                      horzcat(contaminantData.origId, contaminantData.origProcId), 'rows', 'stable');
                  
    index = find(~particleData.active{(impactTime - 1),1});
    impactPrevious = intersect(horzcat(particleData.origId{(impactTime - 1),1}(index,1), particleData.origProcId{(impactTime - 1),1}(index,1)), ...
                      horzcat(contaminantData.origId, contaminantData.origProcId), 'rows', 'stable');
                  
    impactsOfInterest = setdiff(impactsCurrent, impactPrevious, 'rows', 'stable');
    
    index = find(ismember(horzcat(particleData.origId{end,1}, particleData.origProcId{end,1}), ...
                 impactsOfInterest, 'rows'));
             
    % Collate Particles of Interest
    clearvars contaminantData impactsCurrent impactsPrevious impactsOfInterest;
    
    for i = 1:height(particleProps)
        prop = particleProps{i,1};
        contaminantData.(prop) = particleData.(prop){end,1}(index,:);
    end
    
end

disp(' ');

% Limit Tracking Count
switch method

    case {'A', 'B'}
        disp([num2str(height(contaminantData.active)), ' Valid Particles Recorded on Surface']);

    case 'C'
        disp([num2str(height(contaminantData.active)), ' Valid Particles Recorded in the Far-Field']);

end

if height(contaminantData.active) > 500
    disp('    WARNING: Tracking a Large Number of Particles Will Be Computationally Expensive and Difficult to Visualise');
    disp('             It Is Recommended to Track No More Than 500 Particles');
end

valid = false;
while ~valid
    disp(' ');
    selection = input('Limit Tracking Count? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        count = height(contaminantData.active);
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        count = inputCount();

        if count > height(contaminantData.active)
            disp('    WARNING: Value Exceeds Available Particle Count');
        else
            valid = true;
        end

    else
        disp('        WARNING: Invalid Entry');
    end

end

% Select a Random Set of Particles to Track
if count ~= height(contaminantData.active)
    index = sort(randperm(height(contaminantData.active), count))';

    for i = 1:height(particleProps)
        prop = particleProps{i,1};
        contaminantData.(prop) = contaminantData.(prop)(index,:);
    end

end

disp(' ');

disp('***********');
disp('  Running  ');

tic;

disp(' ');

% Perform Tracking
disp(['    Tracking ', num2str(height(contaminantData.active)), ' Particles']);

for i = 1:height(particleData.time)    
    trackingData.ID = intersect(horzcat(particleData.origId{i,1}, particleData.origProcId{i,1}), horzcat(contaminantData.origId, contaminantData.origProcId), 'rows', 'stable');
    index = find(ismember(horzcat(particleData.origId{i,1}, particleData.origProcId{i,1}), trackingData.ID, 'rows'));

    for j = 1:height(trackingData.ID)
        trackingData.position{j,i} = particleData.positionCartesian{i,1}(index(j,1),:);
    end

end

for i = 1:height(trackingData.position)
    trackingData.path{i,1} = unique(cell2mat(trackingData.position(i,:)'), 'stable', 'rows');
    
    if interp && height(trackingData.path{i,1}) >= 3
        trackingData.path{i,1} = interparc((round(deltaT / 1e-3) * height(trackingData.path{i,1})), ...
                                           trackingData.path{i,1}(:,1), ...
                                           trackingData.path{i,1}(:,2), ...
                                           trackingData.path{i,1}(:,3), 'spline');
    end
    
end

disp(' ');

% Plot Particle Tracks
disp('    Plotting Particle Tracks');

switch method
    
    case 'A'
        trackingType = 'Base';
        
    case 'B'
        trackingType = 'Body';
        
    case 'C'
        trackingType = 'Far_Field';
        
end

if exist('impactTime', 'var')
    timeInst = erase(num2str(str2double(particleData.time{impactTime,1}), '%.4f'), '.');
    figName = ['Particle_Tracking', trackingType, '_T', timeInst, '_D', num2str(minD), ...
               '_D', num2str(maxD), '_C', num2str(count)];
else
    figName = ['Particle_Tracking', trackingType, '_D', num2str(minD), '_D', num2str(maxD), ...
               '_C', num2str(count)];
end

cMap = viridis(256);
particlePath = trackingData.path;
particleD = contaminantData.d;
figTitle = {' ', ' '};

fig = trackingPlots(fig, figName, cMap, geometry, particlePath, ...
                    particleD, maxD, minD, figTitle, ...
                    xLims, yLims, zLims);

disp(' ');

executionTime = toc;

disp(['    Run Time: ', num2str(executionTime), 's']);
disp(' ');
disp('  Success  ')
disp('***********');

disp(' ');


%% Cleaning

clearvars -except contaminantData TrackingData impactTime count minD maxD ;
disp(' ');


%% Local Functions

function T = inputTime(time)

    valid = false;
    while ~valid
        T = str2double(input('    Impact Time [s]: ', 's'));

        if isnan(T) || length(T) > 1
            disp('        WARNING: Invalid Entry');
        else
            valid = true;
        end

    end
    
    T = find(ismember(time, num2str(T)));

end


function D = inputD(type)

    valid = false;
    while ~valid
        D = str2double(input(['    ', type, ' Diameter of Interest [', char(956), 'm]: '], 's'));

        if isnan(D) || length(D) > 1
            disp('        WARNING: Invalid Entry');
        else
            valid = true;
        end

    end

end


function C = inputCount()

    valid = false;
    while ~valid
        C = str2double(input('    Number of Particles to Track: ', 's'));

        if isnan(C) || length(C) > 1
            disp('        WARNING: Invalid Entry');
        else
            valid = true;
        end

    end

end