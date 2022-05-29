%% Probe Data Reader v1.2
% ----
% Collates and Optionally Saves OpenFOAM v7 Probe Data Generated using 'PODmeshGenerator.m'
% ----
% Usage: data = readProbeData(caseFolder, timeDirs, timePrecision field, nProc);
%        'caseFolder'    -> Case Path Stored as String
%        'caseName'      -> Case Name, Stored as a String
%        'timeDirs'      -> Time Directories, Obtained With 'timeDirectories.m'
%        'timePrecision' -> Required Rounding Precision for 'deltaT', Obtained With 'timeDirectories.m'
%        'field'         -> Desired Field Stored as String
%        'nProc'         -> Number of Processors Used for Parallel Collation


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Formatting Changes
% v1.2 - Minor Update to Support Additional Versatility of 'velocityProcessing.m'


%% Supported OpenFOAM Cases

% Run_Test
% Windsor_2022


%% Supported Fields

% Pressure: 'p'
% Velocity: 'U'


%% Main Function

function data = readProbeData(caseFolder, caseName, timeDirs, timePrecision, field, nProc) %#ok<INUSD>

    if ~isempty(timeDirs)
        disp('Probe Data Identified in the Following Time Directories:');

        for i = 1:height(timeDirs)
            disp(['    ', timeDirs(i).name]);
        end
        
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
                continue
            end
            
            endTime = inputTime('End');
            
            if endTime == -1
                continue
            elseif endTime < startTime
                disp('        WARNING: Invalid Time Format (''endTime'' Precedes ''startTime'')');
                continue
            elseif endTime < str2double(timeDirs(1).name) || startTime > str2double(timeDirs(end).name)
                disp('        WARNING: No Probe Data in Selected Time Range');
                continue
            end

            i = 1;
            while i <= height(timeDirs)
                
                if str2double(timeDirs(i).name) < startTime || str2double(timeDirs(i).name) > endTime
                    timeDirs(i) = [];
                else
                    i = i + 1;
                end
                
            end

            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end

    end
    clear valid;

    % Identify Time Instances
    data.time = zeros(height(timeDirs),1);

    for i = 1:height(timeDirs)
        data.time(i) = (str2double(timeDirs(i).name));
    end

    % Identify Probe Points
    switch field

        case 'p'
            probeType = 'probesPressure';
            
            fileID = fopen([caseFolder, '/postProcessing/', probeType, '/', timeDirs(1).name, '/base_p.xy']);
            data.position = cell2mat(textscan(fileID, '%f %f %f %*[^\n]', 'delimiter', '\n', 'collectOutput', 1));
            [data.position, index] = unique(data.position, 'stable', 'rows'); % Remove Duplicate Entries
            fclose(fileID);

        case 'U'
            probeType = 'probesVelocity';
            
            fileID = fopen([caseFolder, '/postProcessing/', probeType, '/', timeDirs(1).name, '/wake_U.xy']);
            data.position = cell2mat(textscan(fileID, '%f %f %f %*[^\n]', 'delimiter', '\n', 'collectOutput', 1));
            [data.position, index] = unique(data.position, 'stable', 'rows'); % Remove Duplicate Entries
            fclose(fileID);

    end

    disp(' ');

    % Collate Instantaneous Field Data
    disp('***********');
    disp('  Reading  ');

    tic;
    evalc('parpool(nProc);');
    
    % Initialise Progress Bar
    wB = waitbar(0, ['Collating ''', probeType, ''' Data'], 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    dQ = parallel.pool.DataQueue;
    afterEach(dQ, @parforWaitBar);

    parforWaitBar(wB, height(timeDirs));

    switch field

        case 'p'
            p = cell(height(timeDirs),1);

            parfor i = 1:height(timeDirs)
                p{i} = readInstPressureData(caseFolder, probeType, timeDirs, i, index);
                
                send(dQ, []);
            end

            data.p = p;

            clear p;

        case 'U'
            u = cell(height(timeDirs),1);
            v = u;
            w = u;

            parfor i = 1:height(timeDirs)
                [u{i}, v{i}, w{i}] = readInstVelocityData(caseFolder, probeType, timeDirs, i, index);
                
                send(dQ, []);
            end

            data.u = u;
            data.v = v;
            data.w = w;

            clear u v w;

    end
    
    delete(wB);

    evalc('delete(gcp(''nocreate''));');
    executionTime = toc;

    disp(' ');
    
    disp(['    Read Time: ', num2str(executionTime), 's']);

    disp(' ');
    
    % Calculate Time-Averaged Field Data
    disp('    Calculating Time-Averaged Field Data...');

    switch field

        case 'p'
            data.pMean = zeros(height(data.position),1);

            for i = 1:height(timeDirs)
                data.pMean = data.pMean + data.p{i};
            end

            data.pMean = data.pMean / height(timeDirs);

        case 'U'
            data.uMean = zeros(height(data.position),1);
            data.vMean = data.uMean;
            data.wMean = data.uMean;

            for i = 1:height(timeDirs)
                data.uMean = data.uMean + data.u{i};
                data.vMean = data.vMean + data.v{i};
                data.wMean = data.wMean + data.w{i};
            end

            data.uMean = data.uMean / height(timeDirs);
            data.vMean = data.vMean / height(timeDirs);
            data.wMean = data.wMean / height(timeDirs);

    end

    disp(' ');

    % Calculate Instantaneous Fluctuating Field Data
    disp('    Calculating Instantaneous Fluctuating Field Data...');

    switch field

        case 'p'
            data.pPrime = cell(height(timeDirs),1);

            for i = 1:height(timeDirs)
                data.pPrime{i} = data.p{i} - data.pMean;
            end


        case 'U'
            data.uPrime = cell(height(timeDirs),1);
            data.vPrime = data.uPrime;
            data.wPrime = data.uPrime;

            for i = 1:height(timeDirs)
                data.uPrime{i} = data.u{i} - data.uMean;
                data.vPrime{i} = data.v{i} - data.vMean;
                data.wPrime{i} = data.w{i} - data.wMean;
            end

    end
    
    data = orderfields(data);

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
            
            if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/probeData/', caseName, '/', probeType], 'dir')
                mkdir(['/mnt/Processing/Data/Numerical/MATLAB/probeData/', caseName, '/', probeType]);
            end

            startInst = erase(num2str(str2double(timeDirs(1).name), ['%.', num2str(timePrecision), 'f']), '.');
            endInst = erase(num2str(str2double(timeDirs(end).name), ['%.', num2str(timePrecision), 'f']), '.');
            
            fileName = ['/T', startInst, '_T', endInst, '.mat'];
            
            disp(['    Saving to: /mnt/Processing/Data/Numerical/MATLAB/probeData/', caseName, '/', probeType, fileName]);
            save(['/mnt/Processing/Data/Numerical/MATLAB/probeData/', caseName, '/', probeType, fileName], ...
                 'data', '-v7.3', '-noCompression');
            disp('        Success');
            
            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end

    end
    clear valid;

end


%% Local Functions

function time = inputTime(type)

    time = str2double(input(['    Input ', type, ' Time [s]: '], 's'));
    
    if isnan(time) || length(time) > 1 || time <= 0
        disp('        WARNING: Invalid Entry');
        time = -1;
    end

end


function  p = readInstPressureData(caseFolder, probeType, timeDirs, i, index)

    fileID = fopen([caseFolder, '/postProcessing/', probeType, '/', timeDirs(i).name, '/base_p.xy']);
    content = cell2mat(textscan(fileID, '%*f %*f %*f %f', 'delimiter', '\n', 'collectOutput', 1));
    fclose(fileID);

    p = content(index);

end


function  [u, v, w] = readInstVelocityData(caseFolder, probeType, timeDirs, i, index)

    fileID = fopen([caseFolder, '/postProcessing/', probeType, '/', timeDirs(i).name, '/wake_U.xy']);
    content = cell2mat(textscan(fileID, '%*f %*f %*f %f %f %f', 'delimiter', '\n', 'collectOutput', 1));
    fclose(fileID);

    u = content(index,1);
    v = content(index,2);
    w = content(index,3);

end
