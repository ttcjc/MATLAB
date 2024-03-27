%% Surface Lagrangian Data Reader v2.0
% ----
% Collates and Optionally Saves OpenFOAM v7 Surface Lagrangian Data Output
% ----
% Usage: LagData = readLagDataSurface(saveLoc, caseFolder, caseID, dataID, LagProps, ...
%                                     timeDirs, sampleInt, dataFormat);
%
%        'saveLoc'    -> Start of File Path, Stored as a String
%        'caseFolder' -> Case Path, Stored as s String
%        'campaignID' -> Campaign ID, Stored as a String
%        'caseID'     -> Case Name, Stored as a String
%        'dataID'     -> Data ID, Stored as a String
%        'LagProps'   -> Lagrangian Properties to Be Collated, Stored as a Cell Array
%        'timeDirs'   -> Time Directories, Obtained With 'timeDirectories.m'
%        'sampleInt'  -> Data Binning Interval, Must Be a Factor of Original Recording Frequency
%        'dataFormat' -> Data Collation Format, Stored as a String


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Added Support for ‘Snapshot’ Data Collation and the ‘Uslip’ Field
% v2.1 - Sort Particles Based on 'origProcId' and 'origId' to Enable Multiple Injection Sources


%% Suported Data Collation Formats

% Collate All Particle Impacts Between T(n-1) And T(n): 'cumulative'
% Collate Particle Impacts in a Small Window either Side of T(n): 'snapshot'


%% Main Function

function LagData = readLagDataSurface(saveLoc, caseFolder, campaignID, caseID, distributedFiles, dataID, ...
                                      LagProps, timeDirs, sampleInt, dataFormat)
    
    % Collate Planar Lagrangian Data
    disp('============');
    disp('Surface Data');
    disp('============');
    
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

        % Concatenate Distributed Files
        disp('    Concatenating ''LagrangianSurfaceContaminationData'' Files...');

        contentInt = [];
        contentFloat = [];

        for i = 1:height(dataDirs)

            if str2double(dataDirs(i).name) <= str2double(timeDirs(end).name)
                dataFile = dir([caseFolder, '/LagrangianSurfaceContamination/', dataDirs(i).name, '/LagrangianSurfaceContaminationData']);

                fileID = fopen([caseFolder, '/LagrangianSurfaceContamination/', dataDirs(i).name, '/', dataFile.name]);
                contentRaw = textscan(fileID, '%f32 %u32 %u32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32', 'headerLines', 0, 'delimiter', '\n');

                contentInt = [contentInt; cell2mat(contentRaw(2:3))]; %#ok<AGROW>
                contentFloat = [contentFloat; cell2mat(contentRaw([1,(4:end)]))]; %#ok<AGROW>

                fclose(fileID);
            else
                break;
            end

        end
        clear i;

    else

        % Read Data File
        dataFile = dir([caseFolder, '/LagrangianSurfaceContamination/LagrangianSurfaceContaminationData']);

        disp(['    Loading ''', dataFile.name, '''...']);

        fileID = fopen([caseFolder, '/LagrangianSurfaceContamination/', dataFile.name]);
        contentRaw = textscan(fileID, '%f32 %u32 %u32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32', 'headerLines', 0, 'delimiter', '\n');

        contentInt = cell2mat(contentRaw(2:3));
        contentFloat = cell2mat(contentRaw([1,(4:end)]));

        fclose(fileID);
    end
    
    % Reduce Time Instances to Desired Sampling Frequency
    LagData.time = zeros([ceil(height(timeDirs) / sampleInt),1], 'single');
    nTimes = height(LagData.time);

    j = height(timeDirs);
    for i = nTimes:-1:1
        LagData.time(i) = str2double(timeDirs(j).name);
        j = j - sampleInt;
    end
    clear i j;

    % Initialise Particle Properties
    LagData.timeExact = cell(nTimes,1);

    for i = 1:height(LagProps)
        prop = LagProps{i};
        LagData.(prop) = cell(nTimes,1);
    end
    clear i;

    % Collate Particle Data
    disp('        Collating Surface Data...');

    % Initialise Progress Bar
    wB = waitbar(0, 'Collating Surface Data', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';

    % Perform Collation
    switch dataFormat

        case 'cumulative'
            
            timeZero = LagData.time(1) - (LagData.time(2) - LagData.time(1));
            for i = 1:nTimes
                
                if i == 1
                    index = find((contentFloat(:,1) > timeZero) & ...
                                 (contentFloat(:,1) <= LagData.time(1)));
                else
                    index = find((contentFloat(:,1) > LagData.time(i - 1)) & ...
                                 (contentFloat(:,1) <= LagData.time(i)));
                end
                
                if ~isempty(index)
                    LagData.timeExact{i} = contentFloat(index,1);
                    LagData.origId{i} = contentInt(index,1);
                    LagData.origProcId{i} = contentInt(index,2);
                    LagData.d{i} = contentFloat(index,2);
                    LagData.nParticle{i} = contentFloat(index,3);
                    LagData.positionCartesian{i} = contentFloat(index,[4,5,6]);
%                     LagData.U{j} = contentFloat(index,[7,8,9]);
                    LagData.Uslip{i} = contentFloat(index,[10,11,12]);
                end
                
                % Update Waitbar
                waitbar((i / nTimes), wB);
            end
            clear i;

        case 'snapshot'
            
            temporalSmoothing = 100e-6;
            for i = 1:height(LagData.time)
                index = find((contentFloat(:,1) >= (LagData.time(i) - temporalSmoothing)) & ...
                             (contentFloat(:,1) <= LagData.time(i)));
                
                if ~isempty(index)
                    LagData.timeExact{i} = contentFloat(index,1);
                    LagData.origId{i} = contentInt(index,1);
                    LagData.origProcId{i} = contentInt(index,2);
                    LagData.d{i} = contentFloat(index,2);
                    LagData.nParticle{i} = contentFloat(index,3);
                    LagData.positionCartesian{i} = contentFloat(index,[4,5,6]);
%                     LagData.U{j} = contentFloat(index,[7,8,9]);
                    LagData.Uslip{i} = contentFloat(index,[10,11,12]);
                end
                % Update Waitbar
                waitbar((i / nTimes), wB);
            end
            clear i;

    end

    delete(wB);

    disp(' ');
    
    % Sort Particles in ID Order
    disp('    Sorting Particles...');
    
    % Initialise Progress Bar
    wB = waitbar(0, 'Sorting Particles', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    
    % Perform Sort
    for i = 1:nTimes
        [~, index] = sortrows([LagData.origProcId{i}, LagData.origId{i}]);
        
        LagData.timeExact{i} = LagData.timeExact{i}(index);
        
        for j = 1:height(LagProps)
            prop = LagProps{j};
            LagData.(prop){i} = LagData.(prop){i}(index,:);
        end
        clear j;
        
        waitbar((i / nTimes), wB);
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
            
            if ~exist([saveLoc, '/Numerical/MATLAB/LagData/', campaignID, '/', caseID, '/surface'], 'dir')
                mkdir([saveLoc, '/Numerical/MATLAB/LagData/', campaignID, '/', caseID, '/surface']);
            end

            switch dataFormat

                case 'cumulative'
                    disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/LagData/', campaignID, '/', caseID, '/surface/', dataID, '_cumulative.mat']);
                    save([saveLoc, '/Numerical/MATLAB/LagData/', campaignID, '/', caseID, '/surface/', dataID, '_cumulative.mat'], ...
                         'campaignID', 'caseID', 'dataID', 'LagProps', 'LagData', 'sampleInt', 'dataFormat', '-v7.3', '-noCompression');

                case 'snapshot'
                    disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/LagData/', campaignID, '/', caseID, '/surface/', dataID, '_snapshot.mat']);
                    save([saveLoc, '/Numerical/MATLAB/LagData/', campaignID, '/', caseID, '/surface/', dataID, '_snapshot.mat'], ...
                         'campaignID', 'caseID', 'dataID', 'LagProps', 'LagData', 'sampleInt', 'dataFormat', '-v7.3', '-noCompression');

            end
            
            disp('        Success');
            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end
        
    end
    clear valid;
    
end
