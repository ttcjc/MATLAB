%% Time Directory Reader v1.0


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function[timeDirs, deltaT] = timeDirectories(caseFolder)

    disp(['Case: ', caseFolder]);
    disp(' ');
    disp(' ');
    
    if exist([caseFolder, '/0'], 'dir')
        timeDirs = dir(caseFolder);
    else
        error('Invalid Case Directory');
    end

    i = 1;
    while i <= size(timeDirs,1)

        if isnan(str2double(timeDirs(i,1).name))
            timeDirs(i,:) = [];
        else
            i = i + 1;
        end

    end

    deltaT = str2double(timeDirs(end,1).name) - str2double(timeDirs(end-1,1).name);
    
    disp(['Identified ', num2str(size(timeDirs,1)), ' Time Directories']);
    disp(['    ', char(916), 'T = ' num2str(deltaT), 's']);
    
end