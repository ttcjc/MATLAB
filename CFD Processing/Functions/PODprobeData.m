%% POD Probe Data Reader v.1.1
% ----
% Collates and Optionally Saves OpenFOAM v7 Probe Data Generated using 'PODmeshGenerator.m'
% ----
% Usage: [probeData] = PODprobeData2(caseFolder, timeDirs)
%        'caseFolder' -> Case Path Stored as String
%        'timeDirs'   -> Time Directories Identified Using 'timeDirectories.m' 


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Formatting Changes


%% Main Function

function [probeData] = PODprobeData(caseFolder, timeDirs)

    if ~isempty(timeDirs)
        disp('Probe Data Identified in the Following Time Subdirectories:');

        for i = 1:height(timeDirs)
            disp(['    /postProcessing/probesPOD/', timeDirs(i,1).name]);
        end
        
        disp(' ');
    else
        error('No Probe Data Available in Specified Case');
    end
    
    % Select Directories of Interest
    valid = false;
    while ~valid
        disp(' ');
        selection = input('Restrict Data Range? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            valid = true;
            disp(' ');        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            startTime = inputTime('Start');
            endTime = inputTime('End');

            if endTime < str2double(timeDirs(1,1).name) || startTime > str2double(timeDirs(end,1).name)
                disp('        WARNING: No Probe Data in Selected Data Range');
            elseif endTime < startTime
                disp('        WARNING: Invalid Entry');
            else
                i = 1;
                while i <= height(timeDirs)

                    if str2double(timeDirs(i,1).name) < startTime || str2double(timeDirs(i,1).name) > endTime
                        timeDirs(i,:) = [];
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

    for i = 1:size(timeDirs,1)
        probeData.time{i,1} = (str2double(timeDirs(i,1).name));
    end

    probeData.time = cell2mat(probeData.time);

    % Identify Probe Points
    fileID = fopen([caseFolder, '/postProcessing/probesPOD/', timeDirs(1,1).name, '/wake_U.xy']);
    probeData.position = cell2mat(textscan(fileID, '%f %f %f %*[^\n]', 'delimiter', '\n', 'collectOutput', 1));
    [probeData.position, index] = unique(probeData.position, 'stable', 'rows');
    fclose(fileID);

    % Collate Velocity Data
    disp('***********');
    disp('  Reading  ');
    
    tic;
    
    disp(' ');
    
    for i = 1:height(timeDirs)
        disp(['    /postProcessing/probesPOD/', timeDirs(i,1).name, '/wake_U.xy']);
        
        fileID = fopen([caseFolder, '/postProcessing/probesPOD/', timeDirs(i,1).name, '/wake_U.xy']);
        content = cell2mat(textscan(fileID, '%*f %*f %*f %f %f %f', 'delimiter', '\n', 'collectOutput', 1));
        fclose(fileID);

        probeData.u{i,1} = content(index,1);
        probeData.v{i,1} = content(index,2);
        probeData.w{i,1} = content(index,3);
    end
    
    disp(' ');

    % Calculate Time-Average Velocity Components
    disp('    Calculating Time-Average Velocity Components');
    
    probeData.uMean = zeros(size(probeData.position,1),1);
    probeData.vMean = zeros(size(probeData.position,1),1);
    probeData.wMean = zeros(size(probeData.position,1),1);

    for i = 1:height(timeDirs)
        probeData.uMean = probeData.uMean + probeData.u{i,1};
        probeData.vMean = probeData.vMean + probeData.v{i,1};
        probeData.wMean = probeData.wMean + probeData.w{i,1};
    end

    probeData.uMean = probeData.uMean / size(timeDirs,1);
    probeData.vMean = probeData.vMean / size(timeDirs,1);
    probeData.wMean = probeData.wMean / size(timeDirs,1);

    disp(' ');
    
    % Calculate Instantaneous Fluctuating Velocity Components
    disp('    Calculating Instantaneous Fluctuating Velocity Components');
    
    for i = 1:height(timeDirs)
        probeData.uPrime{i,1} = probeData.u{i,1} - probeData.uMean;
        probeData.vPrime{i,1} = probeData.v{i,1} - probeData.vMean;
        probeData.wPrime{i,1} = probeData.w{i,1} - probeData.wMean;
    end
    
    disp(' ');
    
    executionTime = toc;
    
    disp(['    Read Time: ', num2str(executionTime), 's']);
    disp(' ');
    disp('  Success  ')
    disp('***********');
    
    disp(' ');

    valid = false;
    while ~valid
        disp(' ');
        selection = input('Save Probe Data for Future Use? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            namePos = max(strfind(caseFolder, '/')) + 1;
            
            if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/PODprobeData/', caseFolder(namePos(end):end)], 'dir')
                mkdir(['/mnt/Processing/Data/Numerical/MATLAB/PODprobeData/', caseFolder(namePos(end):end)]);
            end
            
            startInst = erase(num2str(str2double(timeDirs(1,1).name), '%.4f'), '.');
            endInst = erase(num2str(str2double(timeDirs(end,1).name), '%.4f'), '.');
            disp(' ');
            disp(['    Saving to: /mnt/Processing/Data/Numerical/MATLAB/PODprobeData/', caseFolder(namePos(end):end), '/T', startInst, '_T', endInst, '.mat']);
            save(['/mnt/Processing/Data/Numerical/MATLAB/PODprobeData/', caseFolder(namePos(end):end), '/T', startInst, '_T', endInst, '.mat'], 'probeData', '-v7.3', '-noCompression');
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