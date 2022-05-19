%% Planar Lagrangian Data Reader v2.0
% ----
% Collates and Optionally Saves OpenFOAM v7 Planar Lagrangian Data Output
% ----
% Usage: LagData = readLagDataPlane(caseFolder, caseName, LagProps, ...
%                                   binInterval, timeDirs, deltaT, timePrecision);
%        'caseFolder'    -> Case Path, Stored as s String
%        'caseName'      -> Case Name, Stored as a String
%        'LagProps'      -> Lagrangian Properties to Be Collated, Stored as a Cell Array
%        'binInterval'   -> Data Binning Interval, Must Be a Factor of Original Recording Frequency
%        'timeDirs'      -> Time Directories, Obtained With 'timeDirectories.m'
%        'deltaT'        -> Time Delta Between Directiories, Obtained With 'timeDirectories.m'
%        'timePrecision' -> Required Rounding Precision for 'deltaT', Obtained With 'timeDirectories.m'


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Parallelised and Split Into Separate Planar, Surface and Volumetric Functions


%% Main Function

function LagData = readLagDataPlane(caseFolder, caseName, LagProps, ...
                                    binInterval, timeDirs, deltaT, timePrecision)
    
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
        LagData.(plane).time = zeros(ceil(height(timeDirs) / binInterval),1);
        
        k = height(timeDirs);
        for j = height(LagData.(plane).time):-1:1
            LagData.(plane).time(j) = str2double(timeDirs(k).name);
            k = k - binInterval;
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
        for j = 1:height(LagData.(plane).time)
            
            if j == 1
                index = find(content(:,1) <= LagData.(plane).time(j));
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
            
%             if ~exist(['~/Data/Numerical/MATLAB/LagData/plane/', caseName], 'dir')
%                 mkdir(['~/Data/Numerical/MATLAB/LagData/plane/', caseName]);
%             end
            
            if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/LagData/plane/', caseName], 'dir')
                mkdir(['/mnt/Processing/Data/Numerical/MATLAB/LagData/plane/', caseName]);
            end
            
            startInst = erase(num2str(str2double(timeDirs(1).name), ['%.', num2str(timePrecision), 'f']), '.');
            endInst = erase(num2str(str2double(timeDirs(end).name), ['%.', num2str(timePrecision), 'f']), '.');
            
            freq = num2str(round((1 / (deltaT * binInterval)), timePrecision));
            
%             save(['~/Data/Numerical/MATLAB/LagData/Volumetric/', caseName, '/T', startInst, '_T', endInst, '_F', freq, '.mat'], 'LagData', 'LagProps', '-v7.3', '-noCompression');
%             disp(['    Saving to: ~/Data/Numerical/MATLAB/LagData/plane/', caseName, '/T', startInst, '_T', endInst, '_F', freq, '.mat']);
%             disp('        Success');

            save(['/mnt/Processing/Data/Numerical/MATLAB/LagData/plane/', caseName, '/T', startInst, '_T', endInst, '_F', freq, '.mat'], 'LagData', 'LagProps', '-v7.3', '-noCompression');
            disp(['    Saving to: /mnt/Processing/Data/Numerical/MATLAB/LagData/plane/', caseName, '/T', startInst, '_T', endInst, '_F', freq, '.mat']);
            disp('        Success');
            
            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end
        
    end
    
end
