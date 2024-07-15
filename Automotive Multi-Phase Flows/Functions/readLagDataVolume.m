%% Volumetric Lagrangian Data Reader v3.1
% ----
% Collates and Optionally Saves OpenFOAM v7 Volumetric Lagrangian Data Output
% ----
% Usage: LagData = readLagDataVolume(saveLoc, caseFolder, caseID, dataID, cloudName, LagProps, ...
%                                    timeDirs, sampleInt, nProc);
%
%        'saveLoc'    -> Start of File Path, Stored as a String
%        'caseFolder' -> Case Path, Stored as s String
%        'campaignID' -> Campaign ID, Stored as a String
%        'caseID'     -> Case Name, Stored as a String
%        'dataID'     -> Data ID, Stored as a String
%        'cloudName'  -> OpenFOAM Cloud Name, Stored as a String
%        'LagProps'   -> Lagrangian Properties to Be Collated, Stored as a Cell Array
%        'timeDirs'   -> Time Directories, Obtained With 'timeDirectories.m'
%        'sampleInt'  -> Data Sample Interval, Must Be a Factor of Original Recording Frequency
%        'nProc'      -> Number of Processors Used for Parallel Collation


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Updated calls to 'globalPos' to 'positionCartesian'
% v2.0 - Rewrite to Support 'LagrangianExtractionPlaneData'
% v3.0 - Parallelised and Split Into Separate Planar, Surface and Volumetric Functions
% v3.1 - Improved Efficiency of Scalar Data Collation and Added Support for the ‘Uslip’ Field
% v3.2 - Sort Particles Based on 'origProcId' and 'origId' to Enable Multiple Injection Sources


%% Main Function

function LagData = readLagDataVolume(saveLoc, caseFolder, campaignID, caseID, dataID, cloudName, LagProps, ...
                                     timeDirs, sampleInt, nProc) %#ok<INUSD>
    
    % Collate Volumetric Lagrangian Data
    disp('===============');
    disp('Volumetric Data');
    disp('===============');
    
    disp(' ');
    
    disp('***********');
    disp(' COLLATING ');
    
    tic;
    
    evalc('parpool(nProc);');
    
    %%%%
    
    disp(' ');
    
    % Reduce Time Instances to Desired Sampling Frequency
    LagData.time = zeros([ceil(height(timeDirs) / sampleInt),1], 'single');
    nTimes = height(LagData.time);

    j = height(timeDirs);
    for i = nTimes:-1:1
        LagData.time(i) = str2double(timeDirs(j).name);
        j = j - sampleInt;
    end
    clear i j;
    
    % Collate Particle Data
    for i = 1:height(LagProps)
        prop = LagProps{i};
        
        disp(['    Collating ''', prop, ''' Data...']);
        
        % Initialise Progress Bar
        wB = waitbar(0, ['Collating ''', prop, ''' Data'], 'name', 'Progress');
        wB.Children.Title.Interpreter = 'none';
        dQ = parallel.pool.DataQueue;
        afterEach(dQ, @parforWaitBar);
        
        parforWaitBar(wB, nTimes);
        
        % Perform Collation
        propData = cell(nTimes,1);
        
        time = LagData.time;
        parfor j = 1:nTimes
            propData{j} = readInstPropData(caseFolder, num2str(time(j), '%.7g'), cloudName, prop);
            
            % Update Waitbar
            send(dQ, []);
        end
        clear time;
        
        delete(wB);
        
        LagData.(prop) = propData; clear propData;
    end
    clear i;

    disp(' ');

    % Sort Particles in ID Order
    disp('    Sorting Particles...');
    
    % Initialise Progress Bar
    wB = waitbar(0, 'Sorting Particles', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    
    % Perform Sort
    for i = 1:nTimes
        [~, index] = sortrows([LagData.origProcId{i}, LagData.origId{i}]);
        
        for j = 1:height(LagProps)
            prop = LagProps{j};
            LagData.(prop){i} = LagData.(prop){i}(index,:);
        end
        clear j;
        
        waitbar((i / nTimes), wB);
    end
    clear i;
    
    delete(wB);
    
    LagData.time = single(LagData.time);
    
    %%%%
    
    evalc('delete(gcp(''nocreate''));');
    
    executionTime = toc;

    disp(' ');

    disp(['    Run Time: ', num2str(executionTime), 's']);

    disp(' ');

    disp('  SUCCESS  ');
    disp('***********');
    
    % Save Data
    valid = false;
    while ~valid
        disp(' ');
        
        selection = input('Save Data for Future Use? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            
            if ~exist([saveLoc, '/Numerical/MATLAB/LagData/', campaignID, '/', caseID, '/volume'], 'dir')
                mkdir([saveLoc, '/Numerical/MATLAB/LagData/', campaignID, '/', caseID, '/volume']);
            end
            
            disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/LagData/', campaignID, '/', caseID, '/volume/', dataID, '.mat']);
            save([saveLoc, '/Numerical/MATLAB/LagData/', campaignID, '/', caseID, '/volume/', dataID, '.mat'], ...
                 'campaignID', 'caseID', 'dataID', 'LagProps', 'LagData', 'sampleInt', '-v7.3', '-noCompression');
            disp('        Success');
            
            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end
        
    end
    clear valid;

end


%% Local Functions

function propData = readInstPropData(caseFolder, time, cloudName, prop)
    
    fileID = fopen([caseFolder, '/', time, '/lagrangian/', cloudName, '/', prop]);
    content = textscan(fileID, '%s', 'headerLines', 7, 'delimiter', '\n');
    
    if contains(content{1}{5}, 'labelField')
        type = 'int';
    else
        type = 'float';
    end        
    
    frewind(fileID);
    
    if height(content{1}) == 13
        
        if contains(content{1}{11}, '{')
            format = 'A';
        elseif contains(content{1}{11}, '((')
            format = 'D';
        else
            format = 'B';
        end
        
    elseif ~isempty(content{1}{11})
        format = 'F';
    elseif contains(content{1}{14}, '(')
        format = 'E';
    else
        format = 'C';
    end
    
    switch format
        
        case 'A'
            dataLine = content{1}{11};
            dataStart = strfind(dataLine, '{') + 1;
            dataEnd = strfind(dataLine, '}') - 1;
            
            nParticles = str2double(dataLine(1:(dataStart - 2)));
            dataVal = str2double(dataLine(dataStart:dataEnd));
            
            switch type
                
                case 'int'
                    propData = uint32(dataVal * ones(nParticles,1));
                    
                case 'float'
                    propData = single(dataVal * ones(nParticles,1));
                    
            end
            
        case 'B'
            dataLine = content{1}{11};
            dataStart = strfind(dataLine, '(') + 1;
            dataEnd = strfind(dataLine, ')') - 1;
            
            dataVals = (dataLine(dataStart:dataEnd));

            switch type

                case 'int'
                    propData = uint64(str2double(split(dataVals, ' ')));
                
                case 'float'
                    propData = single(str2double(split(dataVals, ' ')));
                    
            end
            
        case 'C'
            
            switch type

                case 'int'
                    propData = cell2mat(textscan(fileID, '%u32', 'headerLines', 20, 'delimiter', '\n'));
                
                case 'float'
                    propData = cell2mat(textscan(fileID, '%f32', 'headerLines', 20, 'delimiter', '\n'));
                    
            end
            
        case 'D'
            dataLine = content{1}{11};
            dataStart = strfind(dataLine, '((') + 1;
            dataEnd = strfind(dataLine, '))');
            
            dataVecs = (dataLine(dataStart:dataEnd));
            
            propData = cell2mat(textscan(dataVecs, '(%f32 %f32 %f32)', 'headerLines', 0, 'delimiter', '\n'));
            
        case 'E'
            propData = cell2mat(textscan(fileID, '(%f32 %f32 %f32)', 'headerLines', 20, 'delimiter', '\n'));
        
        case 'F'
            % Barycentric Co-Ordinates Are Confusing and Not Very Useful
            
    end
    
    fclose(fileID);

end
