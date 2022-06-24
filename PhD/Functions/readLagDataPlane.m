%% Planar Lagrangian Data Reader v2.0
% ----
% Collates and Optionally Saves OpenFOAM v7 Planar Lagrangian Data Output
% ----
% Usage: LagData = readLagDataPlane(saveLocation, caseFolder, caseName, dataID, LagProps, ...
%                                   sampleInterval, timeDirs);
%        'saveLocation'  -> Start of File Path, Stored as a String
%        'caseFolder'     -> Case Path, Stored as s String
%        'caseName'       -> Case Name, Stored as a String
%        'dataID'         -> Data ID, Stored as a String
%        'LagProps'       -> Lagrangian Properties to Be Collated, Stored as a Cell Array
%        'sampleInterval' -> Data Binning Interval, Must Be a Factor of Original Recording Frequency
%        'timeDirs'       -> Time Directories, Obtained With 'timeDirectories.m'


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Parallelised and Split Into Separate Planar, Surface and Volumetric Functions


%% Main Function

function LagData = readLagDataPlane(saveLocation, caseFolder, caseName, dataID, LagProps, ...
                                    sampleInterval, timeDirs)
    
    % Collate Planar Lagrangian Data
    disp('===========');
    disp('Planar Data');
    disp('===========');
    
    disp(' ');
    
    disp('***********');
    disp(' COLLATING ');
    
    disp(' ');

    dataFiles = dir([caseFolder, '/LagrangianExtractionPlane/LagrangianExtractionPlaneData_*']);
    
    tic;
    for i = 1:height(dataFiles)
        disp(['    Loading ''', dataFiles(i).name, '''...']);
        
        % Read Data File
        content = readmatrix([caseFolder, '/LagrangianExtractionPlane/', dataFiles(i).name], 'fileType', 'text', 'trailingDelimitersRule', 'ignore');
        
        planePos = strfind(dataFiles(i).name, '_') + 1;
        plane = ['X_', erase(dataFiles(i).name(planePos:end), '.')];
        
        % Align Particles With Plane Position
        content(:,6) = str2double(dataFiles(i).name(planePos:end));
        
        % Bin Particle Impacts Into Desired Frequency Windows
        LagData.(plane).time = zeros(ceil(height(timeDirs) / sampleInterval),1);
        
        k = height(timeDirs);
        for j = height(LagData.(plane).time):-1:1
            LagData.(plane).time(j) = str2double(timeDirs(k).name);
            k = k - sampleInterval;
        end
        clear k;
        
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
            end
            
            waitbar((j / height(LagData.(plane).time)), wB);
        end
        clear timeZero
        
        delete(wB);
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

        for j = 1:height(LagData.(plane).time)
            [LagData.(plane).origId{j}, index] = sort(LagData.(plane).origId{j});
            
            LagData.(plane).timeExact{j} = LagData.(plane).timeExact{j}(index);

            for k = 1:height(LagProps)
                prop = LagProps{k};

                if k == 3
                    % Don't Sort 'origId' Twice
                    continue
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
            
            disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/LagData/', caseName, '/plane/', dataID, '.mat']);
            save([saveLocation, '/Numerical/MATLAB/LagData/', caseName, '/plane/', dataID, '.mat'], ...
                 'dataID', 'LagProps', 'LagData', 'sampleInterval', '-v7.3', '-noCompression');
            disp('        Success');
            
            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end
        
    end
    clear valid;
    
end