%% PIV Data Initialisation v1.0
% ----
% Hi
% ----
% Usage: [] = initialisePIVdata();
%        'Baa' -> Moo

%% Changelog

% v1.0 - Initial Commit


%% Main Function

function [caseFolder, caseName, PIVdata, format, samplingFrequency] = initialisePIVdata(saveLocation, nProc) %#ok<INUSD>

    % Load Previously Collated Data (If Desired/Possible)
    disp('PIV Data Load');
    disp('-----------------------------');
    valid = false;
    while ~valid
        disp(' ');
        selection = input('Load Previously Collated PIV Data? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            break;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2> 
            [fileName, filePath] = uigetfile([saveLocation, '/Experimental/MATLAB/*.mat'], ...
                                             'Select PIV Data');

            if contains(filePath, '/PIV/')
                disp(['    Loading ''', fileName, '''...']);
                PIVdata = load([filePath, fileName], 'PIVdata').PIVdata;
                format = load([filePath, fileName], 'format').format;
                samplingFrequency = load([filePath, fileName], 'samplingFrequency').samplingFrequency;
                disp('        Success');

                caseFolder = [];
                caseName = fileName(1:(end - 4));
                
                return;
            else
                disp('    WARNING: Invalid File Selection');
                clear fileName filePath;
            end

        else
            disp('    WARNING: Invalid Entry');
            clear fileName filePath;
        end

    end

    disp(' ');
    disp(' ');

    % Acquire PIV Data
    disp('PIV Data Acquisition');
    disp('------------------------------------');
    
    caseFolder = uigetdir('~/Data/', 'Select Case');
    
    namePos = max(strfind(caseFolder, '/')) + 1;
    caseName = caseFolder(namePos:end);
    
    disp(' ');

    disp(['Case: ', caseName]);

    % Confirm Support
    if contains(caseFolder, '/Tomographic/')
        format = 'tomo';
    else
        error('Invalid Case Directory (Unsupported Case Type)');
    end
    
    % Confirm Data Availability
    dataFiles = dir([caseFolder, '/*.dat']);
    
    if isempty(dataFiles)
        error('Invalid Case Directory (No PIV Data Found)');
    end

    disp(' ');

    disp('***********');
    disp(' COLLATING ');
    
    tic;
    evalc('parpool(nProc);');
    
    disp(' ');
    
    disp('    Initialising...');
    
    % Specify Sampling Frequency
    samplingFrequency = str2double(caseName((max(strfind(caseName, '_')) + 1):end));
    
    switch format

        case 'planar'
            % Identify Plane Position
            % Initialise Position Grid
            % Convert to Metres and Adjust Origin
            % Sort Position Grid for Compatibility

        case 'stereo'
            % Identify Plane Position
            % Initialise Position Grid
            % Convert to Metres and Adjust Origin
            % Sort Position Grid for Compatibility

        case 'tomo'            
            % Initialise Position Grid
            fileID = fopen([caseFolder, '/', dataFiles(1).name]);
            content = textscan(fileID, '%f %f %f %f %f %f %f %f', 'headerLines', 4, 'delimiter', ' ', 'CollectOutput', true);
            content = content{1};
            fclose(fileID);
        
            % Rotate Co-Ordinate System and Convert to Metres
            PIVdata.positionGrid = zeros(height(content),3);

            PIVdata.positionGrid(:,1) = -content(:,1) / 1000;
            PIVdata.positionGrid(:,2) = -content(:,3) / 1000;
            PIVdata.positionGrid(:,3) = content(:,2) / 1000;
            
            % Sort Position Grid for Compatibility 
            [PIVdata.positionGrid, index] = sortrows(PIVdata.positionGrid,3);

    end
    
    disp(' ');
    
    % Read Instantaneous Data
    disp('    Collating Instantaneous PIV Data...');
    
    % Initialise Progress Bar
    wB = waitbar(0, 'Collating Instantaneous PIV Data', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    dQ = parallel.pool.DataQueue;
    afterEach(dQ, @parforWaitBar);

    parforWaitBar(wB, height(dataFiles));
    
    switch format
    
        case 'planar'
    
        case 'stereo'
    
        case 'tomo'
            time = zeros(height(dataFiles),1);
            u = cell(height(dataFiles),1);
            v = cell(height(dataFiles),1);
            w = cell(height(dataFiles),1);
            corrCoeff = cell(height(dataFiles),1);
            isValid = cell(height(dataFiles),1);

            parfor i = 1:height(dataFiles)
                fileID = fopen([caseFolder, '/', dataFiles(i).name]);
                content = textscan(fileID, '%f %f %f %f %f %f %f %f', 'headerLines', 4, 'delimiter', ' ', 'CollectOutput', true);
                content = content{1};
                fclose(fileID);
        
                time(i) = i * (1 / samplingFrequency);
                u{i} = -content(index,4);
                v{i} = -content(index,6);
                w{i} = content(index,5);
                corrCoeff{i} = content(index,7);
                isValid{i} = content(index,8);
                
                send(dQ, []);
            end

            PIVdata.time = time;
            PIVdata.u = u;
            PIVdata.v = v;
            PIVdata.w = w;
            PIVdata.corrCoeff = corrCoeff;
            PIVdata.isValid = isValid;
            clear time u v w corrCoeff isValid;
    
    end

    delete(wB);

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
            
            saveDir = [saveLocation, '/Experimental/MATLAB/', ...
                       caseFolder((strfind(caseFolder, '/Data/') + 6):(end - length(caseName)))]; 
            
            if ~exist(saveDir, 'dir')
                mkdir(saveDir);
            end
            
            disp(['    Saving to: ', saveDir, caseName, '.mat']);
            save([saveDir, caseName, '.mat'], ...
                 'PIVdata', 'format', 'samplingFrequency', '-v7.3', '-noCompression');
            disp('        Success');
            
            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end
        
    end
    
end