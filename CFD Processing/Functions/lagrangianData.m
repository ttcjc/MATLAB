%% Lagrangian Data Reader v1.0


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function[particleData, particleProps] = lagrangianData(caseFolder, timeDirs)

    % Identify Lagrangian Directories
    i = 1;
    while i <= size(timeDirs,1)

        if ~exist([caseFolder, '/', timeDirs(i,1).name, '/lagrangian/kinematicCloud/active'], 'file')
            timeDirs(i,:) = [];
        else
            i = i + 1;
        end

    end

    if size(timeDirs,1) ~= 0
        disp('Lagrangian Data Identified in the Following Time Subdirectories:');

        for i = 1:size(timeDirs,1)
            disp(['    /', timeDirs(i,1).name]);
        end

        disp(' ');
    else
        error('No Lagrangian Data Available in Specified Case');
    end

    % Select Directories of Interest
    valid = false;
    while ~valid
        disp(' ');
        selection = input('Restrict Data Range? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            valid = true;
            disp(' ');
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            startTime = inputTime('Start');
            endTime = inputTime('End');

            if endTime < str2double(timeDirs(1,1).name) || startTime > str2double(timeDirs(end,1).name)
                disp('        WARNING: No Lagrangian Data in Selected Data Range');
                disp(' ');
            elseif endTime < startTime
                disp('        WARNING: Invalid Entry');
                disp(' ');
            else
                i = 1;
                while i <= size(timeDirs,1)

                    if str2double(timeDirs(i,1).name) < startTime || str2double(timeDirs(i,1).name) > endTime
                        timeDirs(i,:) = [];
                    else
                        i = i + 1;
                    end

                end

                valid = true;
                disp(' ');
            end

        else
            disp('    WARNING: Invalid Entry');
            disp(' ');
        end

    end

    % Select Lagrangian Properties
    particleProps = {'active', 'origId', 'origProcId', 'nParticle', 'positionCartesian', 'd'};

    valid = false;
	while ~valid
        disp(' ');
        selection = input('Store Additional Lagragian Properties? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            valid = true;
            disp(' ');
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            [particlePropsUser, path] = uigetfile([caseFolder, '/', timeDirs(end,1).name, '/lagrangian/kinematicCloud/*.*'], 'Select Additional Lagrangian Properties', 'multiSelect', 'on');

            if contains(path, '/lagrangian/kinematicCloud/')
                valid = true;
                disp(' ');
            else
                clear particlePropsUser path;
                disp('    WARNING: Invalid File Selection');
                disp(' ');
            end

        else
                disp('    WARNING: Invalid Entry');
                disp(' ');
        end

	end
	
	disp(' ');

    disp('Storing the Following Lagarangian Properties:')
	disp('    active            (Required)');
	disp('    origId            (Required)');
	disp('    origProcId        (Required)');
	disp('    nParticle         (Required)');
	disp('    positionCartesian (Required)');
	disp('    d                 (Required)');

    if exist('particlePropsUser', 'var')
        particleProps = [particleProps, particlePropsUser]';

        if isa(particlePropsUser, 'char')
            disp(['    ', particlePropsUser]);
        else

            for i = 1:size(particlePropsUser,2)
                disp(['    ', particlePropsUser{1,i}]);
            end

        end

    else
        particleProps = particleProps';
    end

	for i = 1:size(timeDirs,1)
        particleData.time{i,1} = timeDirs(i,1).name;

        for j = 1:size(particleProps,1)
            prop = particleProps{j,1};
            particleData.(prop) = [];
        end

	end

	disp(' ');
	disp(' ');
	
    % Collect Lagrangian Data
    disp('***********');
    disp('  Reading  ');
    disp(' ');

    tic;
	for i = 1:size(particleData.time,1)
        disp(['    /', particleData.time{i,1}, '/lagrangian/kinematicCloud/']);

        for j = 1:size(particleProps,1)
            prop = particleProps{j,1};
            fileID = fopen([caseFolder, '/' particleData.time{i,1}, '/lagrangian/kinematicCloud/', prop]);
            content = textscan(fileID, '%s', 'headerLines', 15, 'delimiter', '\n', 'collectOutput', 1);

            if size(content{1,1},1) == 5
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
                        particleData.(prop){i,1} = data .* ones(nParticles,3);
					elseif contains(line, '1(')
						dataStart = strfind(line, '(');
						dataEnd = strfind(line, ')');
						nParticles = str2double(line(1:dataStart-1));
						data = str2double(line(dataStart+1:dataEnd-1));
						particleData.(prop){i,1} = data * ones(nParticles,1);
						
                    else
                        dataStart = strfind(line, '{');
                        dataEnd = strfind(line, '}');
                        nParticles = str2double(line(1:dataStart-1));
                        data = str2double(line(dataStart+1:dataEnd-1));
                        particleData.(prop){i,1} = data * ones(nParticles,1);
                    end

                case 'B'
                    fileID = fopen([caseFolder, '/' particleData.time{i,1}, '/lagrangian/kinematicCloud/', prop]);
                    particleData.(prop){i,1} = cell2mat(textscan(fileID, '%f', 'headerLines', 20, 'delimiter', '\n', 'collectOutput', 1));
					fclose(fileID);
					
                case 'C'
                    fileID = fopen([caseFolder, '/' particleData.time{i,1}, '/lagrangian/kinematicCloud/', prop]);              
                    particleData.(prop){i,1} = cell2mat(textscan(fileID, '(%f %f %f %*[^\n]', 'headerLines', 20, 'delimiter', '\n', 'collectOutput', 1));
					fclose(fileID);
					
            end

        end

	end
	
	for i = 1:size(particleData.time,1)
		[particleData.origId{i,1}, index] = sort(particleData.origId{i,1});

		for j = 1:size(particleProps,1)
			prop = particleProps{j,1};

			if j == 2 % Don't re-sort origId
				continue
			else
				particleData.(prop){i,1} = particleData.(prop){i,1}(index,:);
			end

		end
    
	end
    executionTime = toc;

    disp(' ');
    disp(['    Read Time: ', num2str(executionTime), 's']);
    disp(' ');
    disp('  Success  ')
    disp('***********');
    disp(' ');

    valid = false;
    while ~valid
        disp(' ');
        selection = input('Save Particle Data for Future Use? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            namePos = max(strfind(caseFolder, '/'));
			disp(' ');
            disp(['    Saving to: ~/Documents/Engineering/PhD/Data/Numerical/MATLAB/particleData', caseFolder(namePos(end):end), '.mat']);
            save(['~/Documents/Engineering/PhD/Data/Numerical/MATLAB/particleData', caseFolder(namePos(end):end), '.mat'], 'particleData', 'particleProps', '-v7.3', '-noCompression');
            valid = true;
        else
            disp('    WARNING: Invalid Entry');
            disp(' ');
        end

    end

end


%% Local Functions
    
function T = inputTime(type)

    valid = false;
    while ~valid
        disp(' ');
        T = str2double(input(['    ', type, ' Time [s]: '], 's'));

        if isnan(T) || length(T) > 1
            disp('        WARNING: Invalid Entry');
        else
            valid = true;
        end

    end

end
