%% Surface Lagrangian Data Reader v1.0
% ----
% Collates and Optionally Saves OpenFOAM v7 Surface Lagrangian Data Output
% ----
% Usage: LagData = readLagDataSurface(caseFolder, caseName, dataID, LagProps, ...
%                                     sampleInterval, timeDirs);
%        'caseFolder'     -> Case Path, Stored as s String
%        'caseName'       -> Case Name, Stored as a String
%        'dataID'         -> Data ID, Stored as a String
%        'LagProps'       -> Lagrangian Properties to Be Collated, Stored as a Cell Array
%        'sampleInterval' -> Data Binning Interval, Must Be a Factor of Original Recording Frequency
%        'timeDirs'       -> Time Directories, Obtained With 'timeDirectories.m'


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function LagData = readLagDataSurface(caseFolder, caseName, dataID, LagProps, ...
                                      sampleInterval, timeDirs)
    
    % Collate Planar Lagrangian Data
    disp('============');
    disp('Surface Data');
    disp('============');
    
    disp(' ');
    
    disp('***********');
    disp(' COLLATING ');
    
    disp(' ');

    dataFile = dir([caseFolder, '/LagrangianSurfaceContamination/LagrangianSurfaceContaminationData']);
    
    tic;
    
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
    clear j;

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
    for i = 1:height(LagData.time)

        if i == 1
            index = find(content(:,1) <= LagData.time(i));
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
        end

        waitbar((i / height(LagData.time)), wB);
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
                continue
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
            
            if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/LagData/surface/', caseName], 'dir')
                mkdir(['/mnt/Processing/Data/Numerical/MATLAB/LagData/surface/', caseName]);
            end
            
            disp(['    Saving to: /mnt/Processing/Data/Numerical/MATLAB/LagData/surface/', caseName, '/', dataID]);
            save(['/mnt/Processing/Data/Numerical/MATLAB/LagData/surface/', caseName, '/', dataID], ...
                 'dataID', 'LagProps', 'LagData', 'sampleInterval', '-v7.3', '-noCompression');
            disp('        Success');
            
            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end
        
    end
    clear valid;
    
end