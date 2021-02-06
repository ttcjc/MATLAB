%% Lagrangian Data Reader v2.0

% clearvars;
% close all;
% clc;

% caseFolder = ('/home/cjc96/Mount/Uni/OpenFOAM/ttcjc-7/run/Test_Block');
% caseFolder = ('/home/cjc96/Mount/Athena/OpenFOAM/ttcjc-7/run/Windsor_Square_v2');
% [timeDirs, ~] = timeDirectories(caseFolder);
% disp(' ');
% disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Restructured timeDirectories.m and LagrangianData.m Functionality


%% Main Function

function[particleData, nParticles] = lagrangianData(caseFolder, timeDirs)

    % Identify Lagrangian Directories
    i = 1;
    while i <= size(timeDirs,1)
        
        if ~exist([caseFolder, '/', timeDirs(i,1).name, '/lagrangian'], 'dir')
            timeDirs(i,:) = [];
        else
            i = i + 1;
        end
        
    end
    
    if size(timeDirs,1) > 0
        disp('Lagrangian Data Identified in the Following Time Subdirectories:');
        
        for i = 1:size(timeDirs,1)
            disp(['    /', timeDirs(i,1).name]);
        end
        
    else
        error('No Lagrangian Data Available in Specified Directory');
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
            
            i = 1;
            while i <= size(timeDirs,1)
                
                if str2num(timeDirs(i).name) < startTime || str2num(timeDirs(i).name) > endTime %#ok<*ST2NM>
                    timeDirs(i,:) = [];
                else
                    i = i + 1;
                end
                
            end
            
            valid = true;
            disp(' ');
        else
            disp('    WARNING: Invalid Entry');
        end

    end
    
    if size(timeDirs,1) == 0
        error('No Lagrangian Data in Selected Data Range');
    end
    
    % Collect Lagrangian Data
    tic;
    
    particleVars = {'active', 'origId', 'origProcId', 'globalPos', 'd'};
    particleData = cell(2, size(timeDirs,1));
    
    for i = 1:size(timeDirs,1)
        particleData{1,i} = timeDirs(i).name;
        particleData{2,i} = cell(2,size(particleVars,2));
    end
    
    disp('Reading Lagrangian Data:');
    
    dataLength = cell(size(timeDirs,1),1);
    
    for i = 1:size(timeDirs,1)
        disp(['    /', timeDirs(i).name, '/lagrangian/kinematicCloud']);
        
        variableValues = cell(1,size(particleVars,2));
                
        for j = 1:size(particleVars,2)
            fileID = fopen([caseFolder, '/', timeDirs(i).name, '/lagrangian/kinematicCloud/', particleVars{1,j}]);
            dataFile = textscan(fileID, '%s %s %f %f', 'delimiter', ',', 'MultipleDelimsAsOne', 1);
            dataFile = dataFile{1,1};

            if size(dataFile,1) == 18
                dataTMP = dataFile(17,1);
                formatA = strfind(dataTMP{1,1}, '{');
                formatB = strfind(dataTMP{1,1}, '((');

                if isempty(formatA)

                    if isempty(formatB)
                        findOpen = strfind(dataTMP{1,1}, '(');
                        nParticles = str2num(dataTMP{1,1}(1,1:findOpen-1));
                        dataLength{i,1} = nParticles;
                    else
                        findOpen = strfind(dataTMP{1,1}, '(');
                        nParticles = str2num(dataTMP{1,1}(1:findOpen(1)-1));
                        dataLength{i,1} = nParticles;
                    end

                else
                    findOpen = strfind(dataTMP{1,1}, '{');
                    nParticles = str2num(dataTMP{1,1}(1,1:findOpen-1));
                    dataLength{i,1} = nParticles;
                end

            else
                dataTMP = dataFile(19:end-2,1);
                dataLength{i,1} = length(dataTMP);
            end

            fclose(fileID);
        end
        
        for j = 1:size(particleVars,2)
            fileID = fopen([caseFolder, '/', timeDirs(i).name, '/lagrangian/kinematicCloud/', particleVars{1,j}]);
            dataFile = textscan(fileID, '%s %s %f %f', 'delimiter', ',', 'MultipleDelimsAsOne', 1);
            dataFile = dataFile{1,1};
            
            data = cell(dataLength{i,1},1);

            for k = 1:dataLength{i,1}

                if size(dataFile,1) == 18
                    dataTMP = dataFile(17,1);
                    formatA = strfind(dataTMP{1,1}, '{');
                    formatB = strfind(dataTMP{1,1}, '((');
                    
                    if isempty(formatA)

                        if isempty(formatB)
                            findOpen = strfind(dataTMP{1,1}, '(');
                            findClose = strfind(dataTMP{1,1}, ')');
                            nParticles = str2num(dataTMP{1,1}(1,1:findOpen-1));
                            state = str2num(dataTMP{1,1}(findOpen+1:findClose-1));

                            if k <= nParticles
                                data{k,1} = state(k);
                            else
                                break
                            end

                        else
                            findOpen = strfind(dataTMP{1,1}, '(');
                            findClose = strfind(dataTMP{1,1}, ')');
                            nParticles = str2num(dataTMP{1,1}(1,1:findOpen-1));
                            state = str2num(dataTMP{1,1}(1,findOpen(k+1)+1:findClose(k)-1));

                            if k <= nParticles
                                data{k,1} = state;
                            else
                                break
                            end

                        end

                    else
                        findOpen = strfind(dataTMP{1,1}, '{');
                        findClose = strfind(dataTMP{1,1}, '}');
                        nParticles = str2num(dataTMP{1,1}(1,1:findOpen-1));
                        state = str2num(dataTMP{1,1}(findOpen+1:findClose-1));

                        if k <= nParticles
                            data{k,1} = state;
                        else
                            break
                        end

                    end

                elseif dataFile{19,1}(1) == '('
                    dataTMP = dataFile(19:end-2,1);

                    if k <= length(dataTMP)
                        findOpen = strfind(dataTMP{k,1}, '(');
                        findClose = strfind(dataTMP{k,1}, ')');
                        data{k,1} = str2num(dataTMP{k,1}(findOpen+1:findClose-1));
                        data{k,1} = data{k,1}(1:3);
                    else
                        break
                    end

                else
                    dataTMP = dataFile(19:end-2,1);

                    if k <= length(dataTMP)
                        data{k,1} = str2num(dataTMP{k,1});
                    else
                        break
                    end

                end

            end
            
            variableValues{1,j} = data;
            
            particleData{2,i}{1,j} = particleVars{1,j};
            particleData{2,i}{2,j} = variableValues{1,j};
            
            fclose(fileID);      
        end
        
    end

    for i = 1:size(particleData,2)
        dataTMP = cell(length(particleData{2,i}{2,1}),size(particleVars,2));

        for j = 1:size(particleVars,2)
            dataTMP(:,j) = particleData{2,i}{2,j};        
        end

        particleData{2,i} = dataTMP;
        particleData{2,i} = sortrows(particleData{2,i},2);
        particleData{2,i} = [particleVars(1,:) ; particleData{2,i}];
    end

    nParticles = 0;
    for i = 1:size(particleData,2)
        nParticles = max(nParticles,max(cell2mat(particleData{2,i}(2:end,2))));
    end
    nParticles = nParticles+1;

    disp(' ');
    disp(['Identified ', num2str(nParticles), ' Particles']);
    
    executionTime = toc;
    disp(' ');
    disp(['Read Time = ', num2str(executionTime), 's']);
    
    valid = false;
    while ~valid
        disp(' ');
        selection = input('Save Particle Data for Future Use? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            valid = true;
            disp(' ');
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            
            if ~exist('~/MATLAB/CFD_Processing/particleData', 'dir')
                mkdir('~/MATLAB/CFD_Processing/particleData');
            end
    
            namePos = strfind(caseFolder, '/');
            disp(['    Saving to: ~/MATLAB/CFD_Processing/particleData', caseFolder(namePos(end):end), '.mat']);
            save(['~/MATLAB/CFD_Processing/particleData', caseFolder(namePos(end):end), '.mat'], 'particleData');
            valid = true;
            disp(' ');
        else
            disp('    WARNING: Invalid Entry');
        end

    end
    
    disp('*******************************');
    disp('Lagrangian Data Read Successful');
    disp('*******************************');
    
end
    
    
    %% Local Functions
    
    function time = inputTime(type)

        valid = false;
        while ~valid
            time = str2double(input(['    ', type, ' Time [s]: '], 's'));

            if isnan(time)
                disp('    WARNING: Invalid Entry');
                disp(' ');
            elseif length(time) > 1
                disp('    WARNING: Invalid Entry');
                disp(' ');
            else
                valid = true;
            end

        end
    
    end