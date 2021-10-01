%% Planar Lagrangian Data Reader v1.0
% ----
% Collates and Optionally Saves Planar Lagrangian Data Contained In the Custom
% 'LagrangianExtractionPlaneData' File Developed for OpenFOAM v7
% ----
% Usage: [particleDataPlanar, particlePropsPlanar] = lagrangianDataPlanar(caseFolder, timeDirs)
%        'caseFolder' -> Case Path Stored as String
%        'timeDirs'   -> Case Directories as Collated by 'timeDirectories.m'


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function [particleDataPlanar, particlePropsPlanar] = lagrangianDataPlanar(caseFolder, timeDirs)

    % Identify 'LagrangianExtractionPlaneData'
    if ~exist([caseFolder, '/extractionPlaneData/LagrangianExtractionPlaneData'], 'file')
        error('No Lagrangian Data Available in Specified Case');  
    end
    
    % Identify Global Lagrangian Directories
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
    
    disp(' ');
    
    % Collate Planar Lagrangian Data
    disp('***********');
    disp('  Reading  ');
    disp(' ');
    
    tic;

    disp('    Planar Data:');
    disp('        /extractionPlaneData/LagrangianExtractionPlaneData');

    import = readmatrix([caseFolder, '/extractionPlaneData/LagrangianExtractionPlaneData'], ...
                        'numHeaderLines', 9, 'trailingDelimitersRule', 'ignore');

    for i = 1:height(timeDirs)
        particleDataPlanar.time{i,1} = timeDirs(i,1).name;
    end
                    
    for i = 1:height(particleDataPlanar.time)
        
        if i == 1
            index = find(import(:,9) <= str2double(particleDataPlanar.time{i,1}));
        else
            index = find(import(:,9) > str2double(particleDataPlanar.time{i-1,1}) & ...
                         import(:,9) <= str2double(particleDataPlanar.time{i,1}));
        end

        if ~isempty(index)
            particleDataPlanar.positionCartesian{i,1} = import(index,1:3);
            particleDataPlanar.U{i,1} = import(index,4:6);
            particleDataPlanar.d{i,1} = import(index,7);
            particleDataPlanar.nParticle{i,1} = import(index,8);
            particleDataPlanar.impactTime{i,1} = import(index,9);
        end

    end

    particlePropsPlanar = fieldnames(particleDataPlanar);
    particlePropsPlanar(1) = [];
    
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
            
            if ~exist(['~/Documents/Engineering/PhD/Data/Numerical/MATLAB/particleData/', caseFolder(namePos:end)], 'dir')
                mkdir(['~/Documents/Engineering/PhD/Data/Numerical/MATLAB/particleData/', caseFolder(namePos:end)]);
            end
        
            startInst = erase(num2str(str2double(timeDirs(1,1).name), '%.4f'), '.');
            endInst = erase(num2str(str2double(timeDirs(end,1).name), '%.4f'), '.');
            disp(['    Saving to: ~/Documents/Engineering/PhD/Data/Numerical/MATLAB/particleData/', caseFolder(namePos:end), '/T', startInst, '_T', endInst, '_Planar.mat']);
            save(['~/Documents/Engineering/PhD/Data/Numerical/MATLAB/particleData/', caseFolder(namePos:end), '/T', startInst, '_T', endInst, '_Planar.mat'], 'particleDataPlanar', 'particlePropsPlanar', '-v7.3', '-noCompression');
            
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