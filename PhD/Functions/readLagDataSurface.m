%% Surface Lagrangian Data Reader v2.0
% ----
% Collates and Optionally Saves OpenFOAM v7 Surface Lagrangian Data Output
% ----
% Usage: LagData = readLagDataSurface(saveLocation, caseFolder, caseName, dataID, LagProps, ...
%                                     timeDirs, sampleInterval, format);
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
% v2.0 - Added Support for ‘Snapshot’ Data Collation


%% Suported Data Collation Formats

% Collate All Particle Impacts Between T(n-1) And T(n): 'cumulative'
% Collate Particle Impacts in a Small Window either Side of T(n): 'snapshot'


%% Main Function

function LagData = readLagDataSurface(saveLocation, caseFolder, caseName, dataID, LagProps, ...
                                      timeDirs, sampleInterval, format)
    
    % Collate Planar Lagrangian Data
    disp('============');
    disp('Surface Data');
    disp('============');
    
    disp(' ');
    
    disp('***********');
    disp(' COLLATING ');
    
    disp(' ');

    tic;
    
    dataFile = dir([caseFolder, '/LagrangianSurfaceContamination/LagrangianSurfaceContaminationData']);

    disp(['    Loading ''', dataFile.name, '''...']);
    
    % Read Data File
    content = readmatrix([caseFolder, '/LagrangianSurfaceContamination/', dataFile.name], 'fileType', 'text', 'trailingDelimitersRule', 'ignore');

    % Bin Particle Impacts Into Desired Frequency Windows
    LagData.time = zeros(ceil(height(timeDirs) / sampleInterval),1);

    j = height(timeDirs);
    for i = height(LagData.time):-1:1
        LagData.time(i) = str2double(timeDirs(j).name);
        j = j - sampleInterval;
    end

    % Initialise Particle Properties
    LagData.timeExact = cell(height(LagData.time),1);

    for i = 1:height(LagProps)
        prop = LagProps{i};
        LagData.(prop) = cell(height(LagData.time),1);
    end

    % Initialise Progress Bar
    wB = waitbar(0, 'Collating Surface Contamination Data', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';

    % Collate Data
    switch format

        case 'cumulative'
            
            timeZero = LagData.time(1) - (LagData.time(2) - LagData.time(1));
            for i = 1:height(LagData.time)
                
                if i == 1
                    index = find((content(:,1) > timeZero) & ...
                                 (content(:,1) <= LagData.time(1)));
                else
                    index = find((content(:,1) > LagData.time(i - 1)) & ...
                                 (content(:,1) <= LagData.time(i)));
                end
                
                if ~isempty(index)
                    LagData.timeExact{i} = content(index,1);
                    LagData.origId{i} = content(index,2);
                    LagData.origProcId{i} = content(index,3);
                    LagData.d{i} = content(index,4);
                    LagData.nParticle{i} = content(index,5);
                    LagData.positionCartesian{i} = content(index,[6,7,8]);
                    LagData.U{i} = content(index,[9,10,11]);
                    LagData.Uslip{i} = content(index,[12,13,14]);
                end
                
                waitbar((i / height(LagData.time)), wB);
            end

        case 'snapshot'
            
            temporalSmoothing = 5e-5;
            for i = 1:height(LagData.time)
                index = find((content(:,1) >= (LagData.time(i) - temporalSmoothing)) & ...
                             (content(:,1) <= (LagData.time(i) + temporalSmoothing)));
                
                if ~isempty(index)
                    LagData.timeExact{i} = content(index,1);
                    LagData.origId{i} = content(index,2);
                    LagData.origProcId{i} = content(index,3);
                    LagData.d{i} = content(index,4);
                    LagData.nParticle{i} = content(index,5);
                    LagData.positionCartesian{i} = content(index,[6,7,8]);
                    LagData.U{i} = content(index,[9,10,11]);
                    LagData.Uslip{i} = content(index,[12,13,14]);
                end
        
                waitbar((i / height(LagData.time)), wB);
            end

    end

    delete(wB);
    
    executionTime = toc;
    
    disp(' ');

    disp(['    Read Time: ', num2str(executionTime), 's']);

    disp(' ');

    % Sort Particles in ID Order
    disp('    Sorting Particles...');
    
    for i = 1:height(LagData.time)
        [LagData.origId{i}, index] = sort(LagData.origId{i});
        
        LagData.timeExact{i} = LagData.timeExact{i}(index);
        
        for j = 1:height(LagProps)
            prop = LagProps{j};
            
            if j == 3
                % Don't Sort 'origId' Twice
                continue;
            else
                LagData.(prop){i} = LagData.(prop){i}(index,:);
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
            
            if ~exist([saveLocation, '/Numerical/MATLAB/LagData/', caseName, '/surface'], 'dir')
                mkdir([saveLocation, '/Numerical/MATLAB/LagData/', caseName, '/surface']);
            end

            switch format

                case 'cumulative'
                    disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/LagData/', caseName, '/surface/', dataID, '_cumulative.mat']);
                    save([saveLocation, '/Numerical/MATLAB/LagData/', caseName, '/surface/', dataID, '_cumulative.mat'], ...
                    'dataID', 'LagProps', 'LagData', 'sampleInterval', 'format', '-v7.3', '-noCompression');

                case 'snapshot'
                    disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/LagData/', caseName, '/surface/', dataID, '_snapshot.mat']);
                    save([saveLocation, '/Numerical/MATLAB/LagData/', caseName, '/surface/', dataID, '_snapshot.mat'], ...
                    'dataID', 'LagProps', 'LagData', 'sampleInterval', 'format', '-v7.3', '-noCompression');

            end
            
            disp('        Success');
            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end
        
    end
    
end