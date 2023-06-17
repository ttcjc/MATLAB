%% Planar Lagrangian Data Reader v3.0
% ----
% Collates and Optionally Saves OpenFOAM v7 Planar Lagrangian Data Output
% ----
% Usage: LagData = readLagDataPlane(saveLocation, caseFolder, caseID, dataID, LagProps, ...
%                                   timeDirs, sampleInterval, format);
%        'saveLocation'   -> Start of File Path, Stored as a String
%        'caseFolder'     -> Case Path, Stored as s String
%        'caseID'         -> Case Name, Stored as a String
%        'dataID'         -> Data ID, Stored as a String
%        'LagProps'       -> Lagrangian Properties to Be Collated, Stored as a Cell Array
%        'timeDirs'       -> Time Directories, Obtained With 'timeDirectories.m'
%        'sampleInterval' -> Data Binning Interval, Must Be a Factor of Original Recording Frequency
%        'format'         -> Data Collation Format, Stored as a String


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Parallelised and Split Into Separate Planar, Surface and Volumetric Functions
% v3.0 - Added Support for ‘Snapshot’ Data Collation and the ‘Uslip’ Field


%% Suported Data Collation Formats

% Collate All Particle Impacts Between T(n-1) And T(n): 'cumulative'
% Collate Particle Impacts in a Small Window either Side of T(n): 'snapshot'


%% Main Function

function LagData = readLagDataPlane(saveLocation, caseFolder, caseID, distributedFiles, dataID, ...
                                    LagProps, timeDirs, sampleInterval, format)
    
    % Collate Planar Lagrangian Data
    disp('===========');
    disp('Planar Data');
    disp('===========');
    
    disp(' ');
    
    disp('***********');
    disp(' COLLATING ');
    
    tic;
    
    %%%%
    
    disp(' ');
    
    if distributedFiles
        
        % Identify Distributed Directories
        dataDirs = dir([caseFolder, '/LagrangianSurfaceContamination']);

        i = 1;
        while i <= height(dataDirs)
            
            if isnan(str2double(dataDirs(i).name))
                dataDirs(i,:) = [];
            else
                i = i + 1;
            end
            
        end
        clear i;
        
        dataFiles = dir([caseFolder, '/LagrangianExtractionPlane/', dataDirs(1).name, '/LagrangianExtractionPlaneData_*']);
    else
        dataFiles = dir([caseFolder, '/LagrangianExtractionPlane/LagrangianExtractionPlaneData_*']);
    end
    
    % Reduce Time Instances to Desired Sampling Frequency
    sampleTimes = single(zeros(ceil(height(timeDirs) / sampleInterval),1));
    nTimes = height(sampleTimes);

    j = height(timeDirs);
    for i = nTimes:-1:1
        sampleTimes(i) = str2double(timeDirs(j).name);
        j = j - sampleInterval;
    end
    clear i k;
    
    for i = 1:height(dataFiles)
                
        if distributedFiles
            
            % Concatenate Distributed Files
            if i ~= 1
                disp(' ');
            end
            
            disp(['    Concatenating ', dataFiles(i).name, ' Files...']);

            contentInt = [];
            contentFloat = [];

            for j = 1:height(dataDirs)
                
                if str2double(dataDirs(j).name) <= str2double(timeDirs(end).name)
                    fileID = fopen([caseFolder, '/LagrangianExtractionPlane/', dataDirs(j).name, '/', dataFiles(i).name]);
                    contentRaw = textscan(fileID, '%f32 %u32 %u32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32', 'headerLines', 0, 'delimiter', '\n');
                    
                    contentInt = [contentInt; cell2mat(contentRaw(2:3))]; %#ok<AGROW>
                    contentFloat = [contentFloat; cell2mat(contentRaw([1,(4:end)]))]; %#ok<AGROW>
                    
                    fclose(fileID);
                else
                    break;
                end
                
            end
            clear j;
            
        else
            
            % Read Data File
            if i ~= 1
                disp(' ');
            end
            disp(['    Loading ''', dataFiles(i).name, '''...']);
            
            fileID = fopen([caseFolder, '/LagrangianExtractionPlane/', dataFiles(i).name]);
            contentRaw = textscan(fileID, '%f32 %u32 %u32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32', 'headerLines', 0, 'delimiter', '\n');

            contentInt = cell2mat(contentRaw(2:3));
            contentFloat = cell2mat(contentRaw([1,(4:end)]));
            
            fclose(fileID);
        end
        
        % Identify Plane Position
        planePos = strfind(dataFiles(i).name, '_') + 1;

        if strcmp(dataFiles(i).name(planePos), '-')
            plane = ['X_N', erase(dataFiles(i).name((planePos + 1):end), '.')];
        else
            plane = ['X_P', erase(dataFiles(i).name(planePos:end), '.')];
        end
        
        % Align Particles With Plane Position
        contentFloat(:,4) = str2double(dataFiles(i).name(planePos:end));
        
        % Initialise Particle Properties
        LagData.(plane).time = sampleTimes;
        LagData.(plane).timeExact = cell(nTimes,1);
        
        for j = 1:height(LagProps)
            prop = LagProps{j};
            LagData.(plane).(prop) = cell(nTimes,1);
        end
        clear j;
        
        % Collate Particle Data
        disp(['        Collating ''', plane, ''' Data...']);
        
        % Initialise Progress Bar
        wB = waitbar(0, ['Collating ''', plane, ''' Data'], 'name', 'Progress');
        wB.Children.Title.Interpreter = 'none';
        
        % Perform Collation
        switch format

            case 'cumulative'

                timeZero = LagData.(plane).time(1) - (LagData.(plane).time(2) - LagData.(plane).time(1));
                for j = 1:nTimes
                    
                    if j == 1
                        index = find((contentFloat(:,1) > timeZero) & ...
                                     (contentFloat(:,1) <= LagData.(plane).time(1)));
                    else
                        index = find((contentFloat(:,1) > LagData.(plane).time(j - 1)) & ...
                                     (contentFloat(:,1) <= LagData.(plane).time(j)));
                    end
                    
                    if ~isempty(index)
                        LagData.(plane).timeExact{j} = contentFloat(index,1);
                        LagData.(plane).origId{j} = contentInt(index,1);
                        LagData.(plane).origProcId{j} = contentInt(index,2);
                        LagData.(plane).d{j} = contentFloat(index,2);
                        LagData.(plane).nParticle{j} = contentFloat(index,3);
                        LagData.(plane).positionCartesian{j} = contentFloat(index,[4,5,6]);
                        LagData.(plane).U{j} = contentFloat(index,[7,8,9]);
                        LagData.(plane).Uslip{j} = contentFloat(index,[10,11,12]);
                    end
                    
                    waitbar((j / nTimes), wB);
                end
                clear j;

            case 'snapshot'
                
                temporalSmoothing = 100e-6;
                for j = 1:nTimes
                    index = find((contentFloat(:,1) >= (LagData.(plane).time(j) - temporalSmoothing)) & ...
                                 (contentFloat(:,1) <= LagData.(plane).time(j))); 
                
                    if ~isempty(index)
                        LagData.(plane).timeExact{j} = contentFloat(index,1);
                        LagData.(plane).origId{j} = contentInt(index,1);
                        LagData.(plane).origProcId{j} = contentInt(index,2);
                        LagData.(plane).d{j} = contentFloat(index,2);
                        LagData.(plane).nParticle{j} = contentFloat(index,3);
                        LagData.(plane).positionCartesian{j} = contentFloat(index,[4,5,6]);
                        LagData.(plane).U{j} = contentFloat(index,[7,8,9]);
                        LagData.(plane).Uslip{j} = contentFloat(index,[10,11,12]);
                    end
            
                    waitbar((j / nTimes), wB);
                end
                clear j;

        end
        
        delete(wB);
    end
    clear i;

    disp(' ');
    
    % Sort Particles in ID Order
    disp('    Sorting Particles...');

    % Initialise Progress Bar
    wB = waitbar(0, 'Sorting Particles', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    
    % Perform Sort
    for i = 1:height(dataFiles)

        % Identify Plane Position
        planePos = strfind(dataFiles(i).name, '_') + 1;

        if strcmp(dataFiles(i).name(planePos), '-')
            plane = ['X_N', erase(dataFiles(i).name((planePos + 1):end), '.')];
        else
            plane = ['X_P', erase(dataFiles(i).name(planePos:end), '.')];
        end

        for j = 1:nTimes
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
            clear k;
            
            waitbar((j / nTimes), wB);
        end
        clear j;
        
    end
    clear i;
    
    delete(wB);
    
    %%%%

    evalc('delete(gcp(''nocreate''));');
    
    executionTime = toc;

    disp(' ');

    disp(['    Run Time: ', num2str(executionTime), 's']);

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
            
            if ~exist([saveLocation, '/Numerical/MATLAB/LagData/', caseID, '/plane'], 'dir')
                mkdir([saveLocation, '/Numerical/MATLAB/LagData/', caseID, '/plane']);
            end

            switch format

                case 'cumulative'
                    disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/LagData/', caseID, '/plane/', dataID, '_cumulative.mat']);
                    save([saveLocation, '/Numerical/MATLAB/LagData/', caseID, '/plane/', dataID, '_cumulative.mat'], ...
                    'dataID', 'LagProps', 'LagData', 'sampleInterval', 'format', '-v7.3', '-noCompression');

                case 'snapshot'
                    disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/LagData/', caseID, '/plane/', dataID, '_snapshot.mat']);
                    save([saveLocation, '/Numerical/MATLAB/LagData/', caseID, '/plane/', dataID, '_snapshot.mat'], ...
                    'dataID', 'LagProps', 'LagData', 'sampleInterval', 'format', '-v7.3', '-noCompression');

            end
            
            disp('        Success');
            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end
        
    end
    clear valid;
    
end