%% Probe Data Reader v2.0
% ----
% Collates and Optionally Saves OpenFOAM v7 Probe Data Generated using 'PODmeshGenerator.m'
% ----
% Usage: [dataID, probeData, sampleInt] = readProbeData(saveLoc, caseFolder, campaignID, caseID, ...
%                                                       timeDirs, deltaT, timePrecision, ...
%                                                       probeType, probeRegion, nProc);
%
%        'saveLoc'       -> Start of File Path, Stored as a String
%        'caseFolder'    -> Case Path, Stored as s String
%        'campaignID'    -> Campaign ID, Stored as a String
%        'caseID'        -> Case Name, Stored as a String
%        'timeDirs'      -> Time Directories, Obtained With 'timeDirectories.m'
%        'deltaT'        -> Time Delta Between Directiories, Obtained With 'timeDirectories.m'
%        'timePrecision' -> Required Rounding Precision for 'deltaT', Obtained With 'timeDirectories.m'
%        'probeType'     -> Probe Type, Stored as a String (See Supported Probe Types)
%        'probeRegion'   -> Probe Region, Stored as a String (See OpenFOAM case)
%        'nProc'         -> Number of Processors Used for Parallel Collation

%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Formatting Changes
% v1.2 - Minor Update to Support Additional Versatility of 'velocityProcessing.m'
% v2.0 - Update To Improve Consistency of Structures Across Repository


%% Supported OpenFOAM Campaigns

% Windsor_Upstream_2023
% Windsor_fullScale


%% Supported Probe Types

% Pressure: 'pProbes'
% Velocity: 'uProbes'


%% Main Function

%#ok<*INUSD>

function [dataID, probeData, sampleInt] = readProbeData(saveLoc, caseFolder, campaignID, caseID, ...
                                                        timeDirs, deltaT, timePrecision, ...
                                                        probeType, probeRegion, nProc)

    if ~isempty(timeDirs)
        disp('Probe Data Identified in the Following Time Directories:');

        for i = 1:height(timeDirs)
            disp(['    ', timeDirs(i).name]);
        end
        clear i;
        
    else
        error('Invalid Case Directory (No Probe Data Available)');
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
            sampleInt = 1;
            
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            sampleInt = inputFreq(round((1 / deltaT), timePrecision));
            
            if sampleInt == -1
                continue;
            end
            
            if sampleInt >= (floor(height(timeDirs) / 2))
                disp('            WARNING: Sampling Interval Must Fall Within Data Range');
            else
                valid = true;
            end                
            
        else
            disp('        WARNING: Invalid Entry');
        end

    end
    clear valid;
    
    disp(' ');

    % Collate Probe Data
    disp('***********');
    disp(' COLLATING ');

    tic;
    
    evalc('parpool(nProc);');
    
    %%%%
    
    disp(' ');
    
    disp('    Initialising...');
    
    % Define Data ID
    startInst = erase(num2str(str2double(timeDirs(1).name), ['%.', num2str(timePrecision), 'f']), '.');
    endInst = erase(num2str(str2double(timeDirs(end).name), ['%.', num2str(timePrecision), 'f']), '.');
    freq = num2str(round((1 / (deltaT * sampleInt)), timePrecision));

    dataID = ['T', startInst, '_T', endInst, '_F', freq];

    % Reduce Time Instances to Desired Sampling Frequency
    probeData.time = zeros([ceil(height(timeDirs) / sampleInt),1], 'single');
    nTimes = height(probeData.time);
    
    j = height(timeDirs);
    for i = height(probeData.time):-1:1
        probeData.time(i) = str2double(timeDirs(j).name);
        j = j - sampleInt;
    end
    clear i j;

    % Identify Probe Points
    switch probeType

        case 'pProbes'
            fileID = fopen([caseFolder, '/postProcessing/probesPressure/', ...
                           num2str(probeData.time(1), '%.7g'), '/', probeRegion, '_p.xy']);
            
        case 'uProbes'
            fileID = fopen([caseFolder, '/postProcessing/probesVelocity/', ...
                           num2str(probeData.time(1), '%.7g'), '/', probeRegion, '_U.xy']);
            
    end
    
    probeData.positionGrid = cell2mat(textscan(fileID, '%f32 %f32 %f32 %*[^\n]', 'delimiter', '\n'));
    fclose(fileID);
    
    % Remove Duplicate Entries
    [probeData.positionGrid, index] = unique(probeData.positionGrid, 'stable', 'rows');
    
    disp(' ');
    
    % Collate Instantaneous Field Data
    disp('    Collating Instantaneous Field Data...');
    
    % Initialise Progress Bar
    wB = waitbar(0, 'Collating Instantaneous Field Data', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    dQ = parallel.pool.DataQueue;
    afterEach(dQ, @parforWaitBar);

    parforWaitBar(wB, nTimes);
    
    % Perform Collation
    switch probeType

        case 'pProbes'
            p = cell(nTimes,1);

            time = probeData.time;
            parfor i = 1:nTimes
                p{i} = readInstPressureData(caseFolder, num2str(time(i), '%.7g'), ...
                                            probeRegion, index);
                
                % Update Waitbar
                send(dQ, []);
            end
            clear time;
            
            probeData.p.inst = p; clear p;
            
        case 'uProbes'
            u = cell(height(probeData.time),1);
            v = cell(height(probeData.time),1);
            w = cell(height(probeData.time),1);

            time = probeData.time;
            parfor i = 1:nTimes
                [u{i}, v{i}, w{i}] = readInstVelocityData(caseFolder, num2str(time(i), '%.7g'), ...
                                                          probeRegion, index);
                
                % Update Waitbar
                send(dQ, []);
            end
            clear time;
            
            probeData.u.inst = u; clear u;
            probeData.v.inst = v; clear v;
            probeData.w.inst = w; clear w;
            
    end
    
    delete(wB);
    
    % Order Struct Fields
    switch probeType
        
        case 'pProbes'
            probeData = orderfields(probeData, {'positionGrid', 'time', 'p'});
            
        case 'uProbes'
            probeData = orderfields(probeData, {'positionGrid', 'time', 'u', 'v', 'w'});
            
    end
    
    %%%%
    
    evalc('delete(gcp(''nocreate''));');
    
    executionTime = toc;

    disp(' ');
    
    disp(['    Run Time: ', num2str(executionTime), 's']);

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
            
            if ~exist([saveLoc, '/Numerical/MATLAB/probeData/', campaignID, '/', caseID, '/', probeType, '/', probeRegion], 'dir')
                mkdir([saveLoc, '/Numerical/MATLAB/probeData/', campaignID, '/', caseID, '/', probeType, '/', probeRegion]);
            end
            
            disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/probeData/', campaignID, '/', caseID, '/', probeType, '/', probeRegion, '/', dataID, '.mat']);
            
            save([saveLoc, '/Numerical/MATLAB/probeData/', campaignID, '/', caseID, '/', probeType, '/', probeRegion, '/', dataID, '.mat'], ...
                 'campaignID', 'caseID', 'dataID', 'probeType', 'probeRegion', 'probeData', 'sampleInt', 'timePrecision', '-v7.3', '-noCompression');
            
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


function sampleInt = inputFreq(origFreq)
    
    newFreq = str2double(input('        Input Frequency [Hz]: ', 's'));
    
    if isnan(newFreq) || length(newFreq) > 1 || newFreq <= 0
        disp('            WARNING: Invalid Entry');
        
        sampleInt = -1;
    elseif mod(origFreq, newFreq) ~= 0
        disp(['            WARNING: New Frequency Must Be a Factor of ', num2str(origFreq),' Hz']);
        
        sampleInt = -1;
    else
        sampleInt = origFreq / newFreq;
    end
    
end


function  p = readInstPressureData(caseFolder, time, probeRegion, index)

    fileID = fopen([caseFolder, '/postProcessing/probesPressure/', time, '/', probeRegion, '_p.xy']);
    content = cell2mat(textscan(fileID, '%*f32 %*f32 %*f32 %f32', 'delimiter', '\n'));
    fclose(fileID);

    p = content(index);

end


function  [u, v, w] = readInstVelocityData(caseFolder, time, probeRegion, index)

    fileID = fopen([caseFolder, '/postProcessing/probesVelocity/', time, '/', probeRegion, '_U.xy']);
    content = cell2mat(textscan(fileID, '%*f32 %*f32 %*f32 %f32 %f32 %f32', 'delimiter', '\n'));
    fclose(fileID);

    u = content(index,1);
    v = content(index,2);
    w = content(index,3);

end
