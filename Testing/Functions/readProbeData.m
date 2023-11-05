%% Probe Data Reader v1.2
% ----
% Collates and Optionally Saves OpenFOAM v7 Probe Data Generated using 'PODmeshGenerator.m'
% ----
% Usage: [dataID, probeData, sampleInterval] = readProbeData(saveLocation, caseFolder, caseID, timeDirs, ...
%                                                            deltaT, timePrecision, probeType, nProc)
%
%        'saveLocation'  -> Start of File Path, Stored as a String
%        'caseFolder'    -> Case Path, Stored as a String
%        'caseID'        -> Case Name, Stored as a String
%        'timeDirs'      -> Time Directories, Obtained With 'timeDirectories.m'
%        'deltaT'        -> Time Difference Between 'timeDirs', Obtained With 'timeDirectories.m'
%        'timePrecision' -> Required Rounding Precision for 'deltaT', Obtained With 'timeDirectories.m'
%        'probeType'     -> Desired Probe Type, Stored as a String
%        'nProc'         -> Number of Processors Used for Parallel Collation


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Formatting Changes
% v1.2 - Minor Update to Support Additional Versatility of 'velocityProcessing.m'


%% Supported OpenFOAM Cases

% Run_Test
% Windsor_2022


%% Supported Probe Types

% Pressure: 'probesPressure'
% Velocity: 'probesVelocity'


%% Main Function

function [dataID, probeData, sampleInterval] = readProbeData(saveLocation, caseFolder, caseID, timeDirs, ...
                                                             deltaT, timePrecision, probeType, nProc) %#ok<INUSD>

    if ~isempty(timeDirs)
        disp('Probe Data Identified in the Following Time Directories:');

        for i = 1:height(timeDirs)
            disp(['    ', timeDirs(i).name]);
        end
        clear i;
        
    else
        error('No Probe Data Available in Specified Case');
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
                continue;
            end
            
            endTime = inputTime('End');
            
            if endTime == -1
                continue;
            elseif endTime < startTime
                disp('        WARNING: Invalid Time Format (''endTime'' Precedes ''startTime'')');
                continue;
            elseif endTime < str2double(timeDirs(1).name) || startTime > str2double(timeDirs(end).name)
                disp('        WARNING: No Probe Data in Selected Time Range');
                continue;
            end

            i = 1;
            while i <= height(timeDirs)
                
                if str2double(timeDirs(i).name) < startTime || str2double(timeDirs(i).name) > endTime
                    timeDirs(i) = [];
                else
                    i = i + 1;
                end
                
            end
            clear i;

            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end

    end
    clear valid;
    
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
                continue;
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
    clear valid;
    
    % Define Data ID
    startInst = erase(num2str(str2double(timeDirs(1).name), ['%.', num2str(timePrecision), 'f']), '.');
    endInst = erase(num2str(str2double(timeDirs(end).name), ['%.', num2str(timePrecision), 'f']), '.');
    freq = num2str(round((1 / (deltaT * sampleInterval)), timePrecision));

    dataID = ['T', startInst, '_T', endInst, '_F', freq];

    % Reduce Time Instances to Desired Sampling Frequency
    probeData.time = zeros([ceil(height(timeDirs) / sampleInterval),1]);

    j = height(timeDirs);
    for i = height(probeData.time):-1:1
        probeData.time(i) = str2double(timeDirs(j).name);
        j = j - sampleInterval;
    end
    clear i j;

    % Identify Probe Points
    switch probeType

        case 'probesPressure'
            fileID = fopen([caseFolder, '/postProcessing/', probeType, '/', num2str(probeData.time(1), '%.7g'), '/base_p.xy']);
            probeData.position = cell2mat(textscan(fileID, '%f %f %f %*[^\n]', 'delimiter', '\n'));
            [probeData.position, index] = unique(probeData.position, 'stable', 'rows'); % Remove Duplicate Entries
            fclose(fileID);

        case 'probesVelocity'
            fileID = fopen([caseFolder, '/postProcessing/', probeType, '/', num2str(probeData.time(1), '%.7g'), '/wake_U.xy']);
            probeData.position = cell2mat(textscan(fileID, '%f %f %f %*[^\n]', 'delimiter', '\n'));
            [probeData.position, index] = unique(probeData.position, 'stable', 'rows'); % Remove Duplicate Entries
            fclose(fileID);

    end
    
    disp(' ');

    % Collate Instantaneous Field Data
    disp('***********');
    disp('  Reading  ');

    tic;
    
    evalc('parpool(nProc);');
    
    %%%%
    
    nTimes = height(probeData.time);
    
    % Initialise Progress Bar
    wB = waitbar(0, ['Collating ''', probeType, ''' Data'], 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    dQ = parallel.pool.DataQueue;
    afterEach(dQ, @parforWaitBar);

    parforWaitBar(wB, nTimes);

    switch probeType

        case 'probesPressure'
            p = cell(nTimes,1);

            time = probeData.time;
            parfor i = 1:nTimes
                p{i} = readInstPressureData(caseFolder, probeType, num2str(time(i), '%.7g'), index);
                
                % Update Waitbar
                send(dQ, []);
            end
            clear i time;

            probeData.p = p;
            clear p;

        case 'probesVelocity'
            u = cell(height(probeData.time),1);
            v = cell(height(probeData.time),1);
            w = cell(height(probeData.time),1);

            time = probeData.time;
            parfor i = 1:nTimes
                [u{i}, v{i}, w{i}] = readInstVelocityData(caseFolder, probeType, num2str(time(i), '%.7g'), index);
                
                % Update Waitbar
                send(dQ, []);
            end
            clear i time;

            probeData.u = u;
            probeData.v = v;
            probeData.w = w;
            clear u v w;

    end
    
    delete(wB);
    
    %%%%
    
    evalc('delete(gcp(''nocreate''));');
    
    executionTime = toc;

    disp(' ');
    
    disp(['    Read Time: ', num2str(executionTime), 's']);

    disp(' ');
    
    % Calculate Time-Averaged Field Data
    disp('    Calculating Time-Averaged Field Data...');

    switch probeType

        case 'probesPressure'
            probeData.pMean = zeros([height(probeData.position),1]);

            for i = 1:nTimes
                probeData.pMean = probeData.pMean + probeData.p{i};
            end

            probeData.pMean = probeData.pMean / nTimes;

        case 'probesVelocity'
            probeData.uMean = zeros([height(probeData.position),1]);
            probeData.vMean = probeData.uMean;
            probeData.wMean = probeData.uMean;

            for i = 1:nTimes
                probeData.uMean = probeData.uMean + probeData.u{i};
                probeData.vMean = probeData.vMean + probeData.v{i};
                probeData.wMean = probeData.wMean + probeData.w{i};
            end

            probeData.uMean = probeData.uMean / nTimes;
            probeData.vMean = probeData.vMean / nTimes;
            probeData.wMean = probeData.wMean / nTimes;

    end

    disp(' ');

    % Calculate Instantaneous Fluctuating Field Data
    disp('    Calculating Instantaneous Fluctuating Field Data...');

    switch probeType

        case 'probesPressure'
            probeData.pPrime = cell(nTimes,1);

            for i = 1:nTimes
                probeData.pPrime{i} = probeData.p{i} - probeData.pMean;
            end


        case 'probesVelocity'
            probeData.uPrime = cell(nTimes,1);
            probeData.vPrime = probeData.uPrime;
            probeData.wPrime = probeData.uPrime;

            for i = 1:nTimes
                probeData.uPrime{i} = probeData.u{i} - probeData.uMean;
                probeData.vPrime{i} = probeData.v{i} - probeData.vMean;
                probeData.wPrime{i} = probeData.w{i} - probeData.wMean;
            end

    end
    
    probeData = orderfields(probeData);

    disp(' ');

    disp('  Success  ')
    disp('***********');

    % Save Data
    valid = false;
    while ~valid
        disp(' ');
        selection = input('Save Probe Data for Future Use? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            
            if ~exist([saveLocation, '/Numerical/MATLAB/probeData/', caseID, '/', probeType], 'dir')
                mkdir([saveLocation, '/Numerical/MATLAB/probeData/', caseID, '/', probeType]);
            end
            
            disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/probeData/', caseID, '/', probeType, '/', dataID, '.mat']);
            save([saveLocation, '/Numerical/MATLAB/probeData/', caseID, '/', probeType, '/', dataID, '.mat'], ...
                  'caseID', 'dataID', 'probeData', 'sampleInterval', 'timePrecision', '-v7.3', '-noCompression');
            disp('        Success');
            
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


function  p = readInstPressureData(caseFolder, probeType, time, index)

    fileID = fopen([caseFolder, '/postProcessing/', probeType, '/', time, '/base_p.xy']);
    content = cell2mat(textscan(fileID, '%*f %*f %*f %f', 'delimiter', '\n'));
    fclose(fileID);

    p = content(index,4);

end


function  [u, v, w] = readInstVelocityData(caseFolder, probeType, time, index)

    fileID = fopen([caseFolder, '/postProcessing/', probeType, '/', time, '/wake_U.xy']);
    content = cell2mat(textscan(fileID, '%*f %*f %*f %f %f %f', 'delimiter', '\n'));
    fclose(fileID);

    u = content(index,4);
    v = content(index,5);
    w = content(index,6);

end