%% Volumetric Lagrangian Data Reader v1.0
% ----
% Collates and Optionally Saves OpenFOAM v7 Volumetric Lagrangian Data Output
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

function [particleData, particleProps] = readLagrangianDataVolume(caseFolder, timeDirs, deltaT, timePrecision, cloudName, nProc) %#ok<INUSD>
    
    % Confirm Data Availability
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
    
    % Specify Sampling Frequency
    valid = false;
    while ~valid
        disp(' ');
        disp(['Default Sampling Frequency: ', num2str(round((1 / deltaT), timePrecision)), ' Hz']);
        selection = input('    Reduce Recording Frequency? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            sampleInterval = 1;
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            sampleInterval = inputFreq(round((1 / deltaT), timePrecision));
            
            if sampleInterval == -1
                continue
            end
            
            if sampleInterval >= (floor(height(timeDirs) / 2))
                disp('            WARNING: Sampling Interval Must Fall Within Data Range');
            else
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
    
    % Collate Volumetric Lagrangian Data
    disp('***********');
    disp('  Reading  ');

    tic;
    evalc('parpool(nProc);');
    
    % Reduce Time Instances to Desired Sampling Frequency
    particleData.time = zeros(ceil(height(timeDirs) / sampleInterval),1);
    
    j = height(timeDirs);
    for i = height(particleData.time):-1:1
        particleData.time(i) = str2double(timeDirs(j).name);
        j = j - sampleInterval;
    end
    
    % Read Particle Properties
    for i = 1:height(particleProps)
        prop = particleProps{i};
        propData = cell(height(particleData.time),1);
        
        % Initialise Progress Bar
        wB = waitbar(0, ['Collating ''', prop, ''' Data...'], 'name', 'Progress');
        dQ = parallel.pool.DataQueue;
        afterEach(dQ, @parforWaitBar);
        
        parforWaitBar(wB, height(particleData.time));      
        
        % Collate Data
        time = particleData.time;
        parfor j = 1:height(particleData.time)
            propData{j} = readInstPropData(caseFolder, num2str(time(j), ['%.', num2str(timePrecision), 'g']), cloudName, prop);
            send(dQ, []);
        end
        
        particleData.(prop) = propData;
        
        delete(wB);
    end
    
    evalc('delete(gcp(''nocreate''));');
    executionTime = toc;

    disp(' ');

    disp(['    Read Time: ', num2str(executionTime), 's']);

    disp(' ');

    % Sort Particles in ID Order
    disp('    Sorting Particles...');
    
    for i = 1:1:height(particleData.time)
        [particleData.origId{i}, index] = sort(particleData.origId{i});
        
        for j = 1:height(particleProps)
            prop = particleProps{j};
            
            if j == 3
                % Don't Sort 'origId' Twice
                continue
            else
                particleData.(prop){i} = particleData.(prop){i}(index,:);
            end
            
        end
        
    end

    disp(' ');
    
    disp('  Success  ')
    disp('***********');
    
    % Save Data
    valid = false;
    while ~valid
        disp(' ');
        selection = input('Save Particle Data for Future Use? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            namePos = max(strfind(caseFolder, '/')) + 1;

            if ~exist(['~/Data/Numerical/MATLAB/particleData/Volumetric/', caseFolder(namePos(end):end)], 'dir')
                mkdir(['~/Data/Numerical/MATLAB/particleData/Volumetric/', caseFolder(namePos(end):end)]);
            end

            if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/particleData/Volumetric/', caseFolder(namePos(end):end)], 'dir')
                mkdir(['/mnt/Processing/Data/Numerical/MATLAB/particleData/Volumetric/', caseFolder(namePos(end):end)]);
            end

            startInst = erase(num2str(str2double(timeDirs(1).name), ['%.', num2str(timePrecision), 'f']), '.');
            endInst = erase(num2str(str2double(timeDirs(end).name), ['%.', num2str(timePrecision), 'f']), '.');
            
            freq = num2str(round((1 / (deltaT * sampleInterval)), timePrecision));

            save(['~/Data/Numerical/MATLAB/particleData/Volumetric/', caseFolder(namePos(end):end), '/T', startInst, '_T', endInst, '_F', freq, '.mat'], 'particleData', 'particleProps', '-v7.3', '-noCompression');
            disp(['    Saved to: ~/Data/Numerical/MATLAB/particleData/Volumetric/', caseFolder(namePos(end):end), '/T', startInst, '_T', endInst, '_F', freq, '.mat']);

            save(['/mnt/Processing/Data/Numerical/MATLAB/particleData/Volumetric/', caseFolder(namePos(end):end), '/T', startInst, '_T', endInst, '_F', freq, '.mat'], 'particleData', 'particleProps', '-v7.3', '-noCompression');
            disp(['    Saved to: /mnt/Processing/Data/Numerical/MATLAB/particleData/Volumetric/', caseFolder(namePos(end):end), '/T', startInst, '_T', endInst, '_F', freq, '.mat']);

            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end
        
    end

end


%% Local Functions

function time = inputTime(type)

    time = str2double(input(['    Input ', type, ' Time [s]: '], 's'));
    
    if isnan(time) || length(time) > 1 || time <= 0
        disp('        WARNING: Invalid Entry');
        time = -1;
    end

end


function sampleInterval = inputFreq(origFreq)
    
    newFreq = str2double(input('        Input Frequency [Hz]: ', 's'));
    
    if isnan(newFreq) || length(newFreq) > 1 || newFreq <= 0
        disp('            WARNING: Invalid Entry');
        sampleInterval = -1;
    elseif mod(origFreq, newFreq) ~= 0
        disp(['            WARNING: New Frequency Must Be a Factor of ', num2str(origFreq),' Hz']);
        sampleInterval = -1;
    else
        sampleInterval = origFreq / newFreq;
    end
    
end


function propData = readInstPropData(caseFolder, time, cloudName, prop)
    
    fileID = fopen([caseFolder, '/', time, '/lagrangian/', cloudName, '/', prop]);
    content = textscan(fileID, '%s', 'headerLines', 15, 'delimiter', '\n', 'collectOutput', true);
    frewind(fileID);
    
    if height(content{1}) == 5
        
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
            
            propData = dataVal * ones(nParticles,1);
            
        case 'B'
            dataLine = content{1}{3};
            dataStart = strfind(dataLine, '(') + 1;
            dataEnd = strfind(dataLine, ')') - 1;
            
            dataVals = (dataLine(dataStart:dataEnd));
            
            propData = str2double(split(dataVals, ' '));
            
        case 'C'
            dataRange = content{1}(6:(end - 4));
            
            propData = str2double(dataRange);
            
        case 'D'
            dataLine = content{1}{3};
            dataStart = strfind(dataLine, '((') + 1;
            dataEnd = strfind(dataLine, '))');
            
            dataVecs = (dataLine(dataStart:dataEnd));
            
            propData = cell2mat(textscan(dataVecs, '(%f %f %f)', 'delimiter', ' ', 'collectOutput', true));
            
        case 'E'
            propData = cell2mat(textscan(fileID, '(%f %f %f)', 'headerLines', 20, 'delimiter', '\n', 'collectOutput', true));
        
        case 'F'
            % Barycentric Co-Ordinates Are Confusing and Not Very Useful
    
    end
    
    fclose(fileID);

end
    
