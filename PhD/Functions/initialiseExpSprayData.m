%% Experimental Spray Data Initialisation v1.1
% ----
% Initialisation of Experimental Spray Data for Further Processing
% ----
% Usage: [] = [caseFolder, caseID, ...
%              expSprayData, samplingFrequency] = initialiseExpSprayData(saveLocation, nProc)
%
%        'saveLocation'  -> Start of File Path, Stored as a String
%        'nProc'         -> Number of Processors Used for Parallel Collation


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Update to Support Changes to 'plotPlanarScalarField'


%% Main Function

function [caseFolder, caseID, expSprayData, samplingFrequency] = initialiseExpSprayData(saveLocation, nProc) %#ok<INUSD>

    % Load Previously Collated Data (If Desired/Possible)
    disp('Experimental Spray Data Load');
    disp('-----------------------------');
    
    valid = false;
    while ~valid
        disp(' ');
        selection = input('Load Previously Collated Spray Data? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            break;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2> 
            [fileName, filePath] = uigetfile([saveLocation, '/Experimental/MATLAB/planarExperimentalSpray/*.mat'], ...
                                             'Select Experimental Data');

            if contains(filePath, '/planarExperimentalSpray/')
                disp(['    Loading ''', fileName, '''...']);
                
                caseFolder = load([filePath, fileName], 'caseFolder').caseFolder;
                caseID = load([filePath, fileName], 'caseID').caseID;
                expSprayData = load([filePath, fileName], 'expSprayData').expSprayData;
                samplingFrequency = load([filePath, fileName], 'samplingFrequency').samplingFrequency;
                
                disp('        Success');
                
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
    clear valid;

    disp(' ');
    disp(' ');

    % Acquire Experimental Spray Data
    disp('Experimental Spray Data Acquisition');
    disp('------------------------------------');
    
    caseFolder = uigetdir('~/Data', 'Select Case');
    
    namePos = max(strfind(caseFolder, '/')) + 1;
    caseID = caseFolder(namePos:end);
    
    disp(' ');

    disp(['Case: ', caseID]);

    % Confirm Support
    if ~contains(caseFolder, 'Far_Field_Soiling_07_22')
        error('Invalid Case Directory (Unsupported Case Type)');
    end
    
    % Confirm Data Availability
    dataFiles = dir([caseFolder, '/*.xyz']);
    
    if isempty(dataFiles)
        error('Invalid Case Directory (No Spray Data Found)');
    end

    disp(' ');

    disp('***********');
    disp(' COLLATING ');
    
    tic;
    
    evalc('parpool(nProc);');
    
    %%%%
    
    disp(' ');
    
    disp('    Initialising...');
    
    % Specify Sampling Frequency
    samplingFrequency = str2double(caseID((max(strfind(caseID, '_')) + 1):(end - 2)));
    
    % Identify Plane Position
    planePos = str2double(caseID((min(strfind(caseID, '_')) + 1):(strfind(caseID, 'L_') - 1)));
    
    % Initialise Position Grid
    fileID = fopen([caseFolder, '/', dataFiles(1).name]);
    content = cell2mat(textscan(fileID, '%f %f %f', 'headerLines', 0, 'delimiter', '\n'));
    fclose(fileID);

    % Convert to Metres and Adjust Origin
    expSprayData.positionGrid = zeros([height(content),3]);
    expSprayData.positionGrid(:,1) = 0.48325 + (planePos * 1.044);
    expSprayData.positionGrid(:,(2:3)) = content(:,(1:2)) / 1000;
    expSprayData.positionGrid(:,3) = expSprayData.positionGrid(:,3) + 0.1545;
    
    % Sort Position Grid for Compatibility
    [expSprayData.positionGrid, index] = sortrows(expSprayData.positionGrid,3);
    
    disp(' ');
    
    % Read Instantaneous Data
    disp('    Collating Instantaneous Spray Data...');
    
    % Initialise Progress Bar
    wB = waitbar(0, 'Collating Instantaneous Spray Data', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    dQ = parallel.pool.DataQueue;
    afterEach(dQ, @parforWaitBar);

    parforWaitBar(wB, height(dataFiles));
    
    time = zeros([height(dataFiles),1]);
    seedingDensity = cell(height(dataFiles),1);

    parfor i = 1:height(dataFiles)
        fileID = fopen([caseFolder, '/', dataFiles(i).name]);
        content = cell2mat(textscan(fileID, '%f %f %f', 'headerLines', 0, 'delimiter', '\n'));
        fclose(fileID);

        time(i) = i * (1 / samplingFrequency);
        seedingDensity{i} = content(index,3);
        
        % Update Waitbar
        send(dQ, []);
    end
    clear i;
    
    expSprayData.time = time;
    expSprayData.seedingDensity = seedingDensity;
    clear time seedingDensity;
    
    delete(wB);
    
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
            
            if ~exist([saveLocation, '/Experimental/MATLAB/planarExperimentalSpray'], 'dir')
                mkdir([saveLocation, '/Experimental/MATLAB/planarExperimentalSpray']);
            end
            
            disp(['    Saving to: ', saveLocation, '/Experimental/MATLAB/planarExperimentalSpray/', caseID, '.mat']);
            save([saveLocation, '/Experimental/MATLAB/planarExperimentalSpray/', caseID, '.mat'], ...
                 'caseFolder', 'caseID', 'expSprayData', 'samplingFrequency', '-v7.3', '-noCompression');
            disp('        Success');
            
            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end
        
    end
    
end