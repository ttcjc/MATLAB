%% Global Lagrangian Data Reader v2.0
% ----
% Collates and Optionally Saves OpenFOAM v7 Lagrangian Data Output
% ----
% Usage: [particleDataGlobal, particlePropsGlobal] = lagrangianDataGlobal(caseFolder, timeDirs)
%        'caseFolder' -> Case Path Stored as String
%        'timeDirs'   -> Case Directories as Collated by 'timeDirectories.m'


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Updated calls to 'globalPos' to 'positionCartesian'
% v2.0 - Rewrite to Support 'LagrangianExtractionPlaneData'


%% Main Function

function [particleDataGlobal, particlePropsGlobal] = lagrangianDataGlobal(caseFolder, timeDirs)

    % Identify Lagrangian Directories
    i = 1;
    while i <= height(timeDirs)

        if ~exist([caseFolder, '/', timeDirs(i,1).name, '/lagrangian/kinematicCloud/active'], 'file')
            timeDirs(i,:) = [];
        else
            i = i + 1;
        end

    end

    if ~isempty(timeDirs)
        disp('Lagrangian Data Identified in the Following Time Subdirectories:');

        for i = 1:height(timeDirs)
            disp(['    /', timeDirs(i,1).name]);
        end

    else
        error('No Lagrangian Data Available in Specified Case');
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

            if endTime < str2double(timeDirs(1,1).name) || startTime > str2double(timeDirs(end,1).name)
                disp('        WARNING: No Lagrangian Data in Selected Data Range');
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
    
    % Select Lagrangian Properties
    particlePropsGlobal = {'active'; 'origId'; 'origProcId'; 'nParticle'; 'positionCartesian'; 'd'};

    valid = false;
    while ~valid
        disp(' ');
        selection = input('Store Additional Lagragian Properties? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            [particlePropsUser, path] = uigetfile([caseFolder, '/', timeDirs(end,1).name, '/lagrangian/kinematicCloud/*.*'], 'Select Additional Lagrangian Properties', 'multiSelect', 'on');

            if contains(path, '/lagrangian/kinematicCloud/')
                valid = true;
            else
                clear particlePropsUser path;
                disp('    WARNING: Invalid File Selection');
            end

        else
                disp('    WARNING: Invalid Entry');
        end

    end

    disp(' ');
    disp('Storing the Following Lagrangian Properties:')
    disp('    active                (Required)');
    disp('    origId                (Required)');
    disp('    origProcId            (Required)');
    disp('    nParticle             (Required)');
    disp('    positionCartesian     (Required)');
    disp('    d                     (Required)');
    
    if exist('particlePropsUser', 'var')
        particlePropsGlobal = [particlePropsGlobal; particlePropsUser'];

        if isa(particlePropsUser, 'char')
            disp(['    ', particlePropsUser]);
        else

            for i = 1:width(particlePropsUser)
                disp(['    ', particlePropsUser{1,i}]);
            end

        end
        
    end

    for i = 1:height(timeDirs)
        particleDataGlobal.time{i,1} = timeDirs(i,1).name;

        for j = 1:height(particlePropsGlobal)
            prop = particlePropsGlobal{j,1};
            particleDataGlobal.(prop) = [];
        end

    end
    
    disp(' ');
    
    % Collate Global Lagrangian Data
    disp('***********');
    disp('  Reading  ');
    disp(' ');
    
    tic;
    
    disp('    Global Data:');
    
    for i = 1:height(particleDataGlobal.time)
        disp(['        /', particleDataGlobal.time{i,1}, '/lagrangian/kinematicCloud/']);

        for j = 1:height(particlePropsGlobal)
            prop = particlePropsGlobal{j,1};
            fileID = fopen([caseFolder, '/' particleDataGlobal.time{i,1}, '/lagrangian/kinematicCloud/', prop]);
            content = textscan(fileID, '%s', 'headerLines', 15, 'delimiter', '\n', 'collectOutput', 1);

            if height(content{1,1}) == 5
                format = 'A';
            else
                line = textscan(content{1,1}{6,1}, '(%f %f %f', 'delimiter', ' ');

                if isempty(line{1,2})
                    format = 'B';
                else
                    format = 'C';
                end

            end

            fclose(fileID);

            switch format

                case 'A'
                    line = content{1,1}{3,1};

                    if contains(line, '{(')
                        dataStart = strfind(line, '(');
                        dataEnd = strfind(line, ')');
                        nParticles = str2double(line(1:dataStart-2));
                        data = str2double(line(dataStart+1:dataEnd-1));
                        particleDataGlobal.(prop){i,1} = data .* ones(nParticles,3);
                    elseif contains(line, '1(')
                        dataStart = strfind(line, '(');
                        dataEnd = strfind(line, ')');
                        nParticles = str2double(line(1:dataStart-1));
                        data = str2double(line(dataStart+1:dataEnd-1));
                        particleDataGlobal.(prop){i,1} = data * ones(nParticles,1);

                    else
                        dataStart = strfind(line, '{');
                        dataEnd = strfind(line, '}');
                        nParticles = str2double(line(1:dataStart-1));
                        data = str2double(line(dataStart+1:dataEnd-1));
                        particleDataGlobal.(prop){i,1} = data * ones(nParticles,1);
                    end

                case 'B'
                    fileID = fopen([caseFolder, '/' particleDataGlobal.time{i,1}, '/lagrangian/kinematicCloud/', prop]);
                    particleDataGlobal.(prop){i,1} = cell2mat(textscan(fileID, '%f', 'headerLines', 20, 'delimiter', '\n', 'collectOutput', 1));
                    fclose(fileID);

                case 'C'
                    fileID = fopen([caseFolder, '/' particleDataGlobal.time{i,1}, '/lagrangian/kinematicCloud/', prop]);
                    particleDataGlobal.(prop){i,1} = cell2mat(textscan(fileID, '(%f %f %f %*[^\n]', 'headerLines', 20, 'delimiter', '\n', 'collectOutput', 1));
                    fclose(fileID);

            end

        end

    end

    for i = 1:height(particleDataGlobal.time)
        [particleDataGlobal.origId{i,1}, index] = sort(particleDataGlobal.origId{i,1});

        for j = 1:height(particlePropsGlobal)
            prop = particlePropsGlobal{j,1};

            if j == 2 % Don't Sort 'origId' Twice
                continue
            else
                particleDataGlobal.(prop){i,1} = particleDataGlobal.(prop){i,1}(index,:);
            end

        end

    end
    
    executionTime = toc;

    disp(' ');
    disp(['    Read Time: ', num2str(executionTime), 's']);
    disp(' ');
    disp('  Success  ')
    disp('***********');
    
    % Save Lagrangian Data
    valid = false;
    while ~valid
        disp(' ');
        selection = input('Save Particle Data for Future Use? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            namePos = max(strfind(caseFolder, '/')) + 1;
            
            if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/particleData/', caseFolder(namePos:end)], 'dir')
                mkdir(['/mnt/Processing/Data/Numerical/MATLAB/particleData/', caseFolder(namePos:end)]);
            end
        
            startInst = erase(num2str(str2double(timeDirs(1,1).name), '%.4f'), '.');
            endInst = erase(num2str(str2double(timeDirs(end,1).name), '%.4f'), '.');
            disp(['    Saving to: /mnt/Processing/Data/Numerical/MATLAB/particleData/', caseFolder(namePos:end), '/T', startInst, '_T', endInst, '_Global.mat']);
            save(['/mnt/Processing/Data/Numerical/MATLAB/particleData/', caseFolder(namePos:end), '/T', startInst, '_T', endInst, '_Global.mat'], 'particleDataGlobal', 'particlePropsGlobal', '-v7.3', '-noCompression');
            
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
