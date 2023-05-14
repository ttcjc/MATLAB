%% Volumetric Lagrangian Data Reader v3.1
% ----
% Collates and Optionally Saves OpenFOAM v7 Volumetric Lagrangian Data Output
% ----
% Usage: LagData = readLagDataVolume(saveLocation, caseFolder, caseName, dataID, cloudName, LagProps, ...
%                                    timeDirs, sampleInterval, nProc);
%        'saveLocation'   -> Start of File Path, Stored as a String
%        'caseFolder'     -> Case Path, Stored as s String
%        'caseName'       -> Case Name, Stored as a String
%        'dataID'         -> Data ID, Stored as a String
%        'cloudName'      -> OpenFOAM Cloud Name, Stored as a String
%        'LagProps'       -> Lagrangian Properties to Be Collated, Stored as a Cell Array
%        'timeDirs'       -> Time Directories, Obtained With 'timeDirectories.m'
%        'sampleInterval' -> Data Sample Interval, Must Be a Factor of Original Recording Frequency
%        'nProc'          -> Number of Processors Used for Parallel Collation


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Updated calls to 'globalPos' to 'positionCartesian'
% v2.0 - Rewrite to Support 'LagrangianExtractionPlaneData'
% v3.0 - Parallelised and Split Into Separate Planar, Surface and Volumetric Functions
% v3.1 - Improved Efficiency of Scalar Data Collation and Added Support for the ‘Uslip’ Field


%% Main Function

function LagData = readLagDataVolume(saveLocation, caseFolder, caseName, dataID, cloudName, LagProps, ...
                                     timeDirs, sampleInterval, nProc) %#ok<INUSD>
    
    % Collate Volumetric Lagrangian Data
    disp('===============');
    disp('Volumetric Data');
    disp('===============');
    
    disp(' ');
    
    disp('***********');
    disp(' COLLATING ');
    
    tic;
    evalc('parpool(nProc);');
    
    % Reduce Time Instances to Desired Sampling Frequency
    LagData.time = single(zeros(ceil(height(timeDirs) / sampleInterval),1));
    nTimes = height(LagData.time);

    j = height(timeDirs);
    for i = nTimes:-1:1
        LagData.time(i) = str2double(timeDirs(j).name);
        j = j - sampleInterval;
    end
    clear i k;
    
    % Read Particle Properties
    for i = 1:height(LagProps)
        prop = LagProps{i};
        
        disp(['    Collating ''', prop, ''' Data...']);
        
        % Initialise Progress Bar
        wB = waitbar(0, ['Collating ''', prop, ''' Data'], 'name', 'Progress');
        wB.Children.Title.Interpreter = 'none';
        dQ = parallel.pool.DataQueue;
        afterEach(dQ, @parforWaitBar);
        
        parforWaitBar(wB, nTimes);
        
        % Collate Data
        propData = cell(nTimes,1);
        
        time = LagData.time;
        parfor j = 1:nTimes
            propData{j} = readInstPropData(caseFolder, num2str(time(j), '%.7g'), cloudName, prop);
            
            send(dQ, []);
        end
        clear j time;
        
        delete(wB);
        
        LagData.(prop) = propData; clear propData;
    end
    clear i;

    disp(' ');

    % Sort Particles in ID Order
    disp('    Sorting Particles...');
    
    % Initialise Progress Bar
    wB = waitbar(0, 'Sorting Particles', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    
    for i = 1:1:nTimes
        [LagData.origId{i}, index] = sort(LagData.origId{i});
        
        for j = 1:height(LagProps)
            prop = LagProps{j};
            
            if j == 3
                % Don't Sort 'origId' Twice
                continue;
            else
                LagData.(prop){i} = LagData.(prop){i}(index,:);
            end
            
        end
        clear j;
        
        waitbar((i / nTimes), wB);
    end
    clear i;
    
    delete(wB);

    evalc('delete(gcp(''nocreate''));');
    executionTime = toc;

    disp(' ');

    disp(['    Run Time: ', num2str(executionTime), 's']);

    disp(' ');

    disp('  SUCCESS  ');
    disp('***********');
    
    % Save Data
    valid = false;
    while ~valid
        disp(' ');
        selection = input('Save Data for Future Use? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            
            if ~exist([saveLocation, '/Numerical/MATLAB/LagData/', caseName, '/volume'], 'dir')
                mkdir([saveLocation, '/Numerical/MATLAB/LagData/', caseName, '/volume']);
            end
            
            disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/LagData/', caseName, '/volume/', dataID, '.mat']);
            save([saveLocation, '/Numerical/MATLAB/LagData/', caseName, '/volume/', dataID, '.mat'], ...
                 'dataID', 'LagProps', 'LagData', 'sampleInterval', '-v7.3', '-noCompression');
            disp('        Success');
            
            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end
        
    end
    clear valid;

end


%% Local Functions

function propData = readInstPropData(caseFolder, time, cloudName, prop)
    
    fileID = fopen([caseFolder, '/', time, '/lagrangian/', cloudName, '/', prop]);
    content = textscan(fileID, '%s', 'headerLines', 7, 'delimiter', '\n');
    
    if contains(content{1}{5}, 'labelField')
        type = 'int';
    else
        type = 'float';
    end        
    
    frewind(fileID);
    
    if height(content{1}) == 13
        
        if contains(content{1}{3}, '{')
            format = 'A';
        elseif contains(content{1}{3}, '((')
            format = 'D';
        else
            format = 'B';
        end
        
    elseif ~isempty(content{1}{3})
        format = 'F';
    elseif contains(content{1}{6}, '(')
        format = 'E';
    else
        format = 'C';
    end
    
    switch format
        
        case 'A'
            dataLine = content{1}{3};
            dataStart = strfind(dataLine, '{') + 1;
            dataEnd = strfind(dataLine, '}') - 1;
            
            nParticles = str2double(dataLine(1:(dataStart - 2)));
            dataVal = str2double(dataLine(dataStart:dataEnd));
            
            switch type
                
                case 'int'
                    propData = uint32(dataVal * ones(nParticles,1));
                    
                case 'float'
                    propData = single(dataVal * ones(nParticles,1));
                    
            end
            
        case 'B'
            dataLine = content{1}{3};
            dataStart = strfind(dataLine, '(') + 1;
            dataEnd = strfind(dataLine, ')') - 1;
            
            dataVals = (dataLine(dataStart:dataEnd));

            switch type

                case 'int'
                    propData = uint64(str2double(split(dataVals, ' ')));
                
                case 'float'
                    propData = single(str2double(split(dataVals, ' ')));
                    
            end
            
        case 'C'
            
            switch type

                case 'int'
                    propData = cell2mat(textscan(fileID, '%u32', 'headerLines', 20, 'delimiter', '\n'));
                
                case 'float'
                    propData = cell2mat(textscan(fileID, '%f32', 'headerLines', 20, 'delimiter', '\n'));
                    
            end
            
        case 'D'
            dataLine = content{1}{3};
            dataStart = strfind(dataLine, '((') + 1;
            dataEnd = strfind(dataLine, '))');
            
            dataVecs = (dataLine(dataStart:dataEnd));
            
            propData = cell2mat(textscan(dataVecs, '(%f32 %f32 %f32)', 'headerLines', 0, 'delimiter', '\n'));
            
        case 'E'
            propData = cell2mat(textscan(fileID, '(%f32 %f32 %f32)', 'headerLines', 20, 'delimiter', '\n'));
        
        case 'F'
            % Barycentric Co-Ordinates Are Confusing and Not Very Useful
            
    end
    
    fclose(fileID);

end