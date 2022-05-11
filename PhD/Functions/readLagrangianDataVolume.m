%% Volumetric Lagrangian Data Reader v1.0
% ----
% Collates and Optionally Saves OpenFOAM v7 Volumetric Lagrangian Data Output
% ----
% Usage: [particleData, particleProps] = readLagrangianDataVolume(caseFolder, timeDirs, cloudName, nProc)
%        'caseFolder' -> Case Path Stored as String
%        'timeDirs'   -> Time Directories Identified Using 'timeDirectories.m'
%        'cloudName'  -> OpenFOAM Cloud Name Stored as String
%        'nProc'      -> Number of Processors Used for Parallel Collation


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function [particleData, particleProps] = readLagrangianDataVolume(caseFolder, timeDirs, cloudName, nProc) %#ok<INUSD>
    
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
            endTime = inputTime('End');

            if endTime < str2double(timeDirs(1).name) || startTime > str2double(timeDirs(end).name)
                disp('        WARNING: No Lagrangian Data in Selected Time Range');
            elseif endTime < startTime
                disp('        WARNING: Invalid Entry');
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
    
    % Select Lagrangian Properties
    particleProps = {'d'; 'nParticle'; 'origId'; 'origProcId'; 'positionCartesian'; 'U'};
    
    % -> Add Support for Additional Properties

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

    particleData.time = zeros(height(timeDirs),1);
    
    for i = 1:height(timeDirs)
        particleData.time(i) = (str2double(timeDirs(i).name));
    end
    
    tic;
    evalc('parpool(nProc);');
    
    for i = 1:height(particleProps)
        prop = particleProps{i};
        propData = cell(height(timeDirs),1);
        
        wB = waitbar(0, ['Collating ''', prop, ''' Data...']);
        dQ = parallel.pool.DataQueue;
        afterEach(dQ, @parforWaitBar);
        
        parforWaitBar(wB, height(timeDirs));        
        
        parfor j = 1:height(timeDirs)
            propData{j} = readInstPropData(caseFolder, timeDirs, cloudName, prop, j);
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
    
    for i = 1:height(timeDirs)
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

            startInst = erase(num2str(str2double(timeDirs(1).name), '%.4f'), '.');
            endInst = erase(num2str(str2double(timeDirs(end).name), '%.4f'), '.');

            save(['~/Data/Numerical/MATLAB/particleData/Volumetric/', caseFolder(namePos(end):end), '/T', startInst, '_T', endInst, '.mat'], 'particleData', 'particleProps', '-v7.3', '-noCompression');
            disp(['    Saved to: ~/Data/Numerical/MATLAB/particleData/Volumetric/', caseFolder(namePos(end):end), '/T', startInst, '_T', endInst, '.mat']);

            save(['/mnt/Processing/Data/Numerical/MATLAB/particleData/Volumetric/', caseFolder(namePos(end):end), '/T', startInst, '_T', endInst, '.mat'], 'particleData', 'particleProps', '-v7.3', '-noCompression');
            disp(['    Saved to: /mnt/Processing/Data/Numerical/MATLAB/particleData/Volumetric/', caseFolder(namePos(end):end), '/T', startInst, '_T', endInst, '.mat']);

            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end
        
    end

end


%% Local Functions

function T = inputTime(type)

    valid = false;
    while ~valid
        T = str2double(input(['    ', type, ' Time [s]: '], 's'));

        if isnan(T) || length(T) > 1
            disp('        WARNING: Invalid Entry');
        else
            valid = true;
        end

    end

end


function propData = readInstPropData(caseFolder, timeDirs, cloudName, prop, j)
    
    fileID = fopen([caseFolder, '/', timeDirs(j).name, '/lagrangian/', cloudName, '/', prop]);
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
    