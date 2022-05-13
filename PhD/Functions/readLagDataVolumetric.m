%% Volumetric Lagrangian Data Reader v1.0
% ----
% Collates and Optionally Saves OpenFOAM v7 Volumetric Lagrangian Data Output
% ----
% Usage: particleData = readLagDataVolumetric(caseFolder, caseName, cloudName, particleProps, ...
%                                             sampleInterval, timeDirs, deltaT, timePrecision, nProc)
%        'caseFolder'     -> Case Path, Stored as s String
%        'caseName'       -> Case Name, Stored as a String
%        'cloudName'      -> OpenFOAM Cloud Name, Stored as a String
%        'particleProps'  -> Lagrangian Properties to Be Collated, Stored as a Cell Array
%        'sampleInterval' -> Data Sample Interval, Must Be a Factor of Original Recording Frequency
%        'timeDirs'       -> Time Directories, Obtained With 'timeDirectories.m'
%        'deltaT'         -> Time Delta Between Directiories, Obtained With 'timeDirectories.m'
%        'timePrecision'  -> Required Rounding Precision for 'deltaT', Obtained With 'timeDirectories.m'
%        'nProc'          -> Number of Processors Used for Parallel Collation


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function particleData = readLagDataVolumetric(caseFolder, caseName, cloudName, particleProps, ...
                                              sampleInterval, timeDirs, deltaT, timePrecision, nProc) %#ok<INUSD>
    
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
            propData{j} = readInstPropData(caseFolder, num2str(time(j), '%.8g'), cloudName, prop);
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

            if ~exist(['~/Data/Numerical/MATLAB/particleData/volumetric/', caseName], 'dir')
                mkdir(['~/Data/Numerical/MATLAB/particleData/volumetric/', caseName]);
            end

            if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/particleData/volumetric/', caseName], 'dir')
                mkdir(['/mnt/Processing/Data/Numerical/MATLAB/particleData/volumetric/', caseName]);
            end

            startInst = erase(num2str(str2double(timeDirs(1).name), ['%.', num2str(timePrecision), 'f']), '.');
            endInst = erase(num2str(str2double(timeDirs(end).name), ['%.', num2str(timePrecision), 'f']), '.');
            
            freq = num2str(round((1 / (deltaT * sampleInterval)), timePrecision));

%             save(['~/Data/Numerical/MATLAB/particleData/Volumetric/', caseName, '/T', startInst, '_T', endInst, '_F', freq, '.mat'], 'particleData', 'particleProps', '-v7.3', '-noCompression');
%             disp(['    Saved to: ~/Data/Numerical/MATLAB/particleData/volumetric/', caseName, '/T', startInst, '_T', endInst, '_F', freq, '.mat']);

            save(['/mnt/Processing/Data/Numerical/MATLAB/particleData/volumetric/', caseName, '/T', startInst, '_T', endInst, '_F', freq, '.mat'], 'particleData', 'particleProps', '-v7.3', '-noCompression');
            disp(['    Saved to: /mnt/Processing/Data/Numerical/MATLAB/particleData/volumetric/', caseName, '/T', startInst, '_T', endInst, '_F', freq, '.mat']);

            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end
        
    end

end


%% Local Functions

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
    
