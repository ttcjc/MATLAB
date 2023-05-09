%% Planar Lagrangian Data Reader v3.0
% ----
% Collates and Optionally Saves OpenFOAM v7 Planar Lagrangian Data Output
% ----
% Usage: LagData = readLagDataPlane(saveLocation, caseFolder, caseName, dataID, LagProps, ...
%                                   timeDirs, sampleInterval, format);
%        'saveLocation'   -> Start of File Path, Stored as a String
%        'caseFolder'     -> Case Path, Stored as s String
%        'caseName'       -> Case Name, Stored as a String
%        'dataID'         -> Data ID, Stored as a String
%        'LagProps'       -> Lagrangian Properties to Be Collated, Stored as a Cell Array
%        'timeDirs'       -> Time Directories, Obtained With 'timeDirectories.m'
%        'sampleInterval' -> Data Binning Interval, Must Be a Factor of Original Recording Frequency
%        'format'         -> Data Collation Format, Stored as a String


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Parallelised and Split Into Separate Planar, Surface and Volumetric Functions
% v3.0 - Added Support for 'Snapshot' Data Collation


%% Suported Data Collation Formats

% Collate All Particle Impacts Between T(n-1) And T(n): 'cumulative'
% Collate Particle Impacts in a Small Window either Side of T(n): 'snapshot'


%% Main Function

function LagData = readLagDataPlane(saveLocation, caseFolder, caseName, dataID, LagProps, ...
                                              timeDirs, sampleInterval, format)
    
    % Collate Planar Lagrangian Data
    disp('===========');
    disp('Planar Data');
    disp('===========');
    
    disp(' ');
    
    disp('***********');
    disp(' COLLATING ');
    
    disp(' ');

    tic;

    dataFiles = dir([caseFolder, '/LagrangianExtractionPlane/LagrangianExtractionPlaneData_*']);
    
    for i = 1:height(dataFiles)
        disp(['    Loading ''', dataFiles(i).name, '''...']);
        
        % Read Data File
        content = readmatrix([caseFolder, '/LagrangianExtractionPlane/', dataFiles(i).name], 'fileType', 'text', 'trailingDelimitersRule', 'ignore');
        
        % Identify Plane Position
        planePos = strfind(dataFiles(i).name, '_') + 1;

        if strcmp(dataFiles(i).name(planePos), '-')
            plane = ['X_N', erase(dataFiles(i).name((planePos + 1):end), '.')];
        else
            plane = ['X_P', erase(dataFiles(i).name(planePos:end), '.')];
        end
        
        % Align Particles With Plane Position
        content(:,6) = str2double(dataFiles(i).name(planePos:end));
        
        % Bin Particle Impacts Into Desired Frequency Windows
        LagData.(plane).time = zeros(ceil(height(timeDirs) / sampleInterval),1);
        
        k = height(timeDirs);
        for j = height(LagData.(plane).time):-1:1
            LagData.(plane).time(j) = str2double(timeDirs(k).name);
            k = k - sampleInterval;
        end
        
        % Initialise Particle Properties
        LagData.(plane).timeExact = cell(height(LagData.(plane).time),1);
        
        for j = 1:height(LagProps)
            prop = LagProps{j};
            LagData.(plane).(prop) = cell(height(LagData.(plane).time),1);
        end
        
        % Initialise Progress Bar
        wB = waitbar(0, ['Collating ''', plane, ''' Data'], 'name', 'Progress');
        wB.Children.Title.Interpreter = 'none';

        % Collate Data
        switch format

            case 'cumulative'

                timeZero = LagData.(plane).time(1) - (LagData.(plane).time(2) - LagData.(plane).time(1));
                for j = 1:height(LagData.(plane).time)
                    
                    if j == 1
                        index = find((content(:,1) > timeZero) & ...
                                     (content(:,1) <= LagData.(plane).time(1)));
                    else
                        index = find((content(:,1) > LagData.(plane).time(j - 1)) & ...
                                     (content(:,1) <= LagData.(plane).time(j)));
                    end
                    
                    if ~isempty(index)
                        LagData.(plane).timeExact{j} = content(index,1);
                        LagData.(plane).origId{j} = content(index,2);
                        LagData.(plane).origProcId{j} = content(index,3);
                        LagData.(plane).d{j} = content(index,4);
                        LagData.(plane).nParticle{j} = content(index,5);
                        LagData.(plane).positionCartesian{j} = content(index,[6,7,8]);
                        LagData.(plane).U{j} = content(index,[9,10,11]);
                        LagData.(plane).Uslip{j} = content(index,[12,13,14]);
                    end
                    
                    waitbar((j / height(LagData.(plane).time)), wB);
                end

            case 'snapshot'
                
                temporalSmoothing = 5e-5;
                for j = 1:height(LagData.(plane).time)
                    index = find((content(:,1) >= (LagData.(plane).time(j) - temporalSmoothing)) & ...
                                 (content(:,1) <= (LagData.(plane).time(j) + temporalSmoothing)));
                    
                    if ~isempty(index)
                        LagData.(plane).timeExact{j} = content(index,1);
                        LagData.(plane).origId{j} = content(index,2);
                        LagData.(plane).origProcId{j} = content(index,3);
                        LagData.(plane).d{j} = content(index,4);
                        LagData.(plane).nParticle{j} = content(index,5);
                        LagData.(plane).positionCartesian{j} = content(index,[6,7,8]);
                        LagData.(plane).U{j} = content(index,[9,10,11]);
                        LagData.(plane).Uslip{j} = content(index,[12,13,14]);
                    end
            
                    waitbar((j / height(LagData.(plane).time)), wB);
                end

        end
        
        delete(wB);
    end

    executionTime = toc;
    
    disp(' ');

    disp(['    Read Time: ', num2str(executionTime), 's']);

    disp(' ');

    % Sort Particles in ID Order
    disp('    Sorting Particles...');
    
    for i = 1:height(dataFiles)

        % Identify Plane Position
        planePos = strfind(dataFiles(i).name, '_') + 1;

        if strcmp(dataFiles(i).name(planePos), '-')
            plane = ['X_N', erase(dataFiles(i).name((planePos + 1):end), '.')];
        else
            plane = ['X_P', erase(dataFiles(i).name(planePos:end), '.')];
        end

        for j = 1:height(LagData.(plane).time)
            [LagData.(plane).origId{j}, index] = sort(LagData.(plane).origId{j});
            
            LagData.(plane).timeExact{j} = LagData.(plane).timeExact{j}(index);

            for k = 1:height(LagProps)
                prop = LagProps{k};

                if k == 3
                    % Don't Sort 'origId' Twice
                    continue;
                else
                    LagData.(plane).(prop){j} = LagData.(plane).(prop){j}(index,:);
                end

            end

        end
        
    end

    disp(' ');
    
    disp('  SUCCESS  ')
    disp('***********');

    % Save Data
    valid = false;
    while ~valid
        disp(' ');
        selection = input('Save Data for Future Use? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            
            if ~exist([saveLocation, '/Numerical/MATLAB/LagData/', caseName, '/plane'], 'dir')
                mkdir([saveLocation, '/Numerical/MATLAB/LagData/', caseName, '/plane']);
            end

            switch format

                case 'cumulative'
                    disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/LagData/', caseName, '/plane/', dataID, '_cumulative.mat']);
                    save([saveLocation, '/Numerical/MATLAB/LagData/', caseName, '/plane/', dataID, '_cumulative.mat'], ...
                    'dataID', 'LagProps', 'LagData', 'sampleInterval', 'format', '-v7.3', '-noCompression');

                case 'snapshot'
                    disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/LagData/', caseName, '/plane/', dataID, '_snapshot.mat']);
                    save([saveLocation, '/Numerical/MATLAB/LagData/', caseName, '/plane/', dataID, '._snapshot.mat'], ...
                    'dataID', 'LagProps', 'LagData', 'sampleInterval', 'format', '-v7.3', '-noCompression');

            end
            
            disp('        Success');
            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end
        
    end
    
end