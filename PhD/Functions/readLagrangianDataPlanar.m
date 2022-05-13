%% Planar Lagrangian Data Reader v1.0
% ----
% Collates and Optionally Saves OpenFOAM v7 Planar Lagrangian Data Output
% ----
% Usage: [particleData, particleProps] = readLagrangianDataVolume(caseFolder, timeDirs, deltaT, timePrecision, cloudName, nProc)
%        'caseFolder'    -> Case Path Stored as String
%        'timeDirs'      -> Time Directories, Obtained With 'timeDirectories.m'
%        'deltaT'        -> Time Delta Between Directiories, Obtained With 'timeDirectories.m'
%        'timePrecision' -> Required Rounding Precision for 'deltaT', Obtained With 'timeDirectories.m'
%        'cloudName'     -> OpenFOAM Cloud Name Stored as String
%        'nProc'         -> Number of Processors Used for Parallel Collation


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function [particleData, particleProps] = readLagrangianDataPlanar(caseFolder, timeDirs, deltaT, timePrecision, cloudName, nProc) %#ok<INUSD>
    
    % Confirm Data Availability
    if ~exist([caseFolder, '/LagrangianExtractionPlane'], 'dir')
        error('Invalid Case Directory (No Planar Lagrangian Data Available)');
    end
    
    dataFiles = dir([caseFolder, '/LagrangianExtractionPlane/LagrangianExtractionPlaneData_*']);
    
    if isempty(dataFiles)
        error('Invalid Case Directory (No Planar Lagrangian Data Available)');
    end
    
    i = 1;
    while i <= height(timeDirs)

        if ~exist([caseFolder, '/', timeDirs(i).name, '/lagrangian/', cloudName, '/active'], 'file')
            timeDirs(i) = [];
        else
            i = i + 1;
        end

    end
    
    if ~isempty(timeDirs)
        disp('Lagrangian Data Identified in the Following Time Directories:');

        for i = 1:height(timeDirs)
            disp(['    ', timeDirs(i).name]);
        end

    else
        error('Invalid Case Directory (No Volumetric Lagrangian Data Available)');
    end
    
    % Select Times of Interest
    valid = false;
    while ~valid
        disp(' ');
        selection = input('Restrict Data Range? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            startTime = inputTime('Start');
            
            if startTime == -1
                continue
            end
            
            endTime = inputTime('End');
            
            if endTime == -1
                continue
            end
            
            if endTime < startTime
                disp('        WARNING: Invalid Time Format (''endTime'' Precedes ''startTime'')');
            elseif endTime < str2double(timeDirs(1).name) || startTime > str2double(timeDirs(end).name)
                disp('        WARNING: No Lagrangian Data in Selected Time Range');
            else
                
                i = 1;
                while i <= height(timeDirs)

                    if str2double(timeDirs(i).name) < startTime || str2double(timeDirs(i).name) > endTime
                        timeDirs(i) = [];
                    else
                        i = i + 1;
                    end

                end
                
                valid = true;
            end

        else
            disp('    WARNING: Invalid Entry');
        end

    end
    
    % Specify Binning Frequency
    valid = false;
    while ~valid
        disp(' ');
        disp(['Default Binning Frequency: ', num2str(round(1 / deltaT, timePrecision)), ' Hz']);
        selection = input('    Reduce Binning Frequency? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            binInterval = 1;
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            binInterval = inputFreq(round(1 / deltaT, timePrecision));
            
            if binInterval ~= -1
                valid = true;
            end
            
        else
            disp('        WARNING: Invalid Entry');
        end

    end
    
    % Select Lagrangian Properties
    particleProps = {'d'; 'nParticle'; 'origId'; 'origProcId'; 'positionCartesian'; 'U'};

    disp(' ');

    disp('Storing the Following Lagrangian Properties:')
    disp('    d                     (Required)');
    disp('    nParticle             (Required)');
    disp('    origId                (Required)');
    disp('    origProcId            (Required)');
    disp('    positionCartesian     (Required)');
    disp('    U                     (Required)');
    
    disp(' ');
    
    % Collate Planar Lagrangian Data
    disp('***********');
    disp('  Reading  ');
    
    disp(' ');
    
    tic;
    evalc('parpool(nProc);');
    for i = 1:height(dataFiles)
        disp(['    Loading ''', dataFiles(i).name, '''...']);
        
        % Read Data File
        content = readmatrix([caseFolder, '/LagrangianExtractionPlane/', dataFiles(i).name], 'fileType', 'text', 'trailingDelimitersRule', 'ignore');
        
        planePos = strfind(dataFiles(i).name, '_') + 1;
        plane = ['X_', erase(dataFiles(i).name(planePos:end), '.')];
        
        % Bin Particle Impacts Into Desired Frequency Windows
        particleData.(plane).time = zeros((height(timeDirs) / binInterval),1);
        
        k = height(timeDirs);
        for j = height(particleData.(plane).time):-1:1
            particleData.(plane).time(j) = str2double(timeDirs(k).name);
            k = k - binInterval;
        end
        
        % Initialise Particle Properties
        origId = cell(height(particleData.(plane).time),1);
        origProcId = origId;
        d = origId;
        nParticle = origId;
        positionCartesian = origId;
        U = origId;
        
        % Initialise Progress Bar
        wB = waitbar(0, ['Collating ''', plane, ''' Data...'], 'name', 'Progress', 'windowStyle', 'docked');
        dQ = parallel.pool.DataQueue;
        afterEach(dQ, @parforWaitBar);

        parforWaitBar(wB, height(particleData.(plane).time));     

        % Collate Data
        parfor j = 1:height(particleData.(plane).time)
            
            if j == 1
                index = find(content(:,1) <= particleData.(plane).time(j)); %#ok<PFBNS>
            else
                index = find((content(:,1) > particleData.(plane).time(j - 1)) & ...
                             (content(:,1) <= particleData.(plane).time(j)));
            end
            
            if ~isempty(index)
                origId{j} = content(index,2);
                origProcId{j} = content(index,3);
                d{j} = content(index,4);
                nParticle{j} = content(index,5);
                positionCartesian{j} = content(index,[6,7,8]);
                U{j} = content(index,[9,10,11]);
            end
            
            send(dQ, []);            
        end
        
        % Store Data
        particleData.(plane).origId = origId;
        particleData.(plane).origProcId = origProcId;
        particleData.(plane).d = d;
        particleData.(plane).nParticle = nParticle;
        particleData.(plane).positionCartesian = positionCartesian;
        particleData.(plane).U = U;
        
        delete(wB);        
    end
    evalc('delete(gcp(''nocreate''));');
    executionTime = toc;
    
    disp(' ');

    disp(['    Read Time: ', num2str(executionTime), 's']);

    disp(' ');

    % Sort Particles in ID Order
    disp('    Sorting Particles...');
    
    for i = 1:height(dataFiles)
        planePos = strfind(dataFiles(i).name, '_') + 1;
        plane = ['X_', erase(dataFiles(i).name(planePos:end), '.')];

        for j = 1:height(particleData.(plane).time)
            [particleData.(plane).origId{j}, index] = sort(particleData.(plane).origId{j});

            for k = 1:height(particleProps)
                prop = particleProps{k};

                if k == 3
                    % Don't Sort 'origId' Twice
                    continue
                else
                    particleData.(plane).(prop){j} = particleData.(plane).(prop){j}(index,:);
                end

            end

        end
        
    end

    disp(' ');
    
    disp('  Success  ')
    disp('***********');
    
end


%% Local Functions

function time = inputTime(type)

    time = str2double(input(['    Input ', type, ' Time [s]: '], 's'));
    
    if isnan(time) || length(time) > 1 || time <= 0
        disp('        WARNING: Invalid Entry');
        time = -1;
    end

end


function binInterval = inputFreq(origFreq)
    
    newFreq = str2double(input('        Input Frequency [Hz]: ', 's'));
    
    if isnan(newFreq) || length(newFreq) > 1 || newFreq <= 0
        disp('            WARNING: Invalid Entry');
        binInterval = -1;
    elseif mod(origFreq, newFreq) ~= 0
        disp(['            WARNING: New Frequency Must Be a Factor of ', num2str(origFreq),' Hz']);
        binInterval = -1;
    else
        binInterval = origFreq / newFreq;
    end
    
end