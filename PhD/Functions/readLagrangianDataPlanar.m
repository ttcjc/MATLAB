%% Planar Lagrangian Data Reader v1.0
% ----
% Collates and Optionally Saves OpenFOAM v7 Planar Lagrangian Data Output
% ----
% Usage: [particleData, particleProps] = readLagrangianDataVolume(caseFolder, timeDirs, cloudName, nProc)
%        'caseFolder' -> Case Path Stored as String
%        'timeDirs'   -> Time Directories Identified Using 'timeDirectories.m'
%        'cloudName'  -> OpenFOAM Cloud Name Stored as String
%        'freq'       -> Desired Binning Frequency [Hz]


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function [particleData, particleProps] = readLagrangianDataPlanar(caseFolder, timeDirs, cloudName, freq)
    
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
            endTime = inputTime('End');

            if endTime < str2double(timeDirs(1).name) || startTime > str2double(timeDirs(end).name)
                disp('        WARNING: No Lagrangian Data in Selected Time Range');
            elseif endTime < startTime
                disp('        WARNING: Invalid Entry');
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
    
    % Select Lagrangian Properties
    particleProps = {'d'; 'nParticle'; 'origId'; 'origProcId'; 'positionCartesian'; 'U'};
    
    % -> Add Support for Additional Properties

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
    for i = 1:height(dataFiles)
        disp(['    ', dataFiles(i).name]);
                
        % Read Data File
        content = readmatrix([caseFolder, '/LagrangianExtractionPlane/', dataFiles(i).name], 'fileType', 'text', 'trailingDelimitersRule', 'ignore');
        
        planePos = strfind(dataFiles(i).name, '_') + 1;
        plane = ['X_', erase(dataFiles(i).name(planePos:end), '.')];
        
        % Bin Particle Impacts Into Desired Frequency Windows
        particleData.(plane).time = (str2double(timeDirs(1).name):(1 / freq):str2double(timeDirs(end).name))';
        
        for j = 1:height(particleProps)
            prop = particleProps{j};
            particleData.(plane).(prop) = cell(height(particleData.(plane).time),1);
        end
        
        for j = 1:height(particleData.(plane).time)
            
            if j == 1
                index = find(content(:,1) <= particleData.(plane).time(j));
            else
                index = find((content(:,1) > particleData.(plane).time(j - 1)) & ...
                             (content(:,1) <= particleData.(plane).time(j)));
            end
            
            if ~isempty(index)
                particleData.(plane).origId{j} = content(index,2);
                particleData.(plane).origProcId{j} = content(index,3);
                particleData.(plane).d{j} = content(index,4);
                particleData.(plane).nParticle{j} = content(index,5);
                particleData.(plane).positionCartesian{j} = content(index,[6,7,8]);
                particleData.(plane).U{j} = content(index,[9,10,11]);
            end
            
        end
        
    end
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

function T = inputTime(type)

    valid = false;
    while ~valid
        T = str2double(input(['    ', type, ' Time [s]: '], 's'));

        if isnan(T) || length(T) > 1
            disp('        WARNING: Invalid Entry');
        else
            valid = true;
        end

    end

end