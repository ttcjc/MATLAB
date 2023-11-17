%% Experimental Spray Data Initialisation v2.0
% ----
% Initialisation of Experimental Spray Data for Further Processing
% ----
% Usage: [campaignID, caseID, planeID, ...
%         expSprayData, samplingFrequency] = initialiseExpSprayData(saveLoc, dataLoc, nProc)
%
%        'saveLoc' -> Start of File Path, Stored as a String
%        'dataLoc' -> Start of File Path, Stored as a String
%        'nProc'   -> Number of Processors Used for Parallel Collation


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Update to Support Changes to 'plotPlanarScalarField'
% v2.0 - Update To Support Changes to DaVis Data Format


%% Main Function

function [campaignID, caseID, planeID, ...
          expSprayData, sampleFreq] = initialiseExpSprayData(saveLoc, dataLoc, nProc) %#ok<INUSD>

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
            [fileName, filePath] = uigetfile([saveLoc, '/Experimental/MATLAB/experimentalSprayData/*.mat'], ...
                                             'Select Experimental Data');

            if ~isnumeric(filePath)

                if contains(filePath, '/experimentalSprayData/')
                    disp(['    Loading ''', fileName, '''...']);

                    campaignID = load([filePath, fileName], 'campaignID').campaignID;
                    caseID = load([filePath, fileName], 'caseID').caseID;
                    planeID = load([filePath, fileName], 'planeID').planeID;
                    expSprayData = load([filePath, fileName], 'expSprayData').expSprayData;
                    sampleFreq = load([filePath, fileName], 'sampleFreq').sampleFreq;

                    disp('        Success');

                    return;
                else
                    disp('    WARNING: Invalid File Selection');
                    clear fileName filePath;
                end

            else
                error('No File Selected');
            end

        else
            disp('    WARNING: Invalid Entry');
        end

    end
    clear valid;

    disp(' ');
    disp(' ');

    % Acquire Experimental Spray Data
    disp('Experimental Spray Data Acquisition');
    disp('------------------------------------');

    valid = false;
    while ~valid
        disp(' ');
        
        disp('Select DaVis Results Folder...');

        caseFolder = uigetdir([dataLoc, '/Experimental'], 'Select DaVis Results Folder');

        if ~isnumeric(caseFolder)

            % Confirm Support
            if contains(caseFolder, 'Far_Field_Soiling_07_22')
                dataFiles = dir([caseFolder, '/*.csv']);
            else
                disp('    WARNING: Invalid Case Directory (Unsupported Case Type)');

                continue;
            end

        else
            error('No Results Folder Selected');
        end

        % Confirm Data Availability
        if isempty(dataFiles)
            disp('    WARNING: Invalid Case Directory (No Spray Data Found)');

            continue;
        else
            valid = true;
        end

    end
    clear valid

    campaignID = caseFolder((strfind(caseFolder, 'Experimental/') + 13):(strfind(caseFolder, '/Results') - 1));
    caseID = caseFolder((max(strfind(caseFolder, '/')) + 1):end);

    disp(' ');

    disp(['Case: ', campaignID, ', ', caseID]);
    
    disp(' ');
    
    disp('***********');
    disp(' COLLATING ');
    
    tic;
    
    %%%%
    
    disp(' ');
    
    disp('    Initialising...');
    
    evalc('parpool(nProc);');
    
    % Identify Sampling Frequency
    sampleFreq = str2double(caseID((strfind(caseID, 's_') + 2):(strfind(caseID, 'Hz') - 1)));
    
    % Identify Plane Position
    planeID = caseID((min(strfind(caseID, '_')) + 1):(strfind(caseID, 'L_')));
    
    % Initialise Position Grid
    fileID = fopen([caseFolder, '/', dataFiles(1).name]);
    content = cell2mat(textscan(fileID, '%f32;%f32;%f32', 'headerLines', 1, 'delimiter', '\n'));
    fclose(fileID);

    % Convert to Metres and Adjust Origin
    expSprayData.positionGrid = zeros([height(content),3], 'single');
    expSprayData.positionGrid(:,1) = 0.48325 + (str2double(planeID(1:(end-1))) * 1.044); % Offset Due to Plane Position
    expSprayData.positionGrid(:,(2:3)) = content(:,(1:2)) / 1000;
    expSprayData.positionGrid(:,3) = expSprayData.positionGrid(:,3) + (0.309 / 2); % Offset Due to 309-15-3 Calibration Board
    
    % Sort Position Grid for 'ndgrid' Compatibility
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
    
    time = zeros([height(dataFiles),1], 'single');
    density = cell(height(dataFiles),1);

    parfor i = 1:height(dataFiles)
        fileID = fopen([caseFolder, '/', dataFiles(i).name]);
        content = cell2mat(textscan(fileID, '%f32;%f32;%f32', 'headerLines', 1, 'delimiter', '\n'));
        fclose(fileID);

        time(i) = i * (1 / sampleFreq);
        density{i} = content(index,3);
        
        % Update Waitbar
        send(dQ, []);
    end
    clear i;
    
    expSprayData.time = time;
    expSprayData.density = density;
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
            
            if ~exist([saveLoc, '/Experimental/MATLAB/experimentalSprayData/', campaignID], 'dir')
                mkdir([saveLoc, '/Experimental/MATLAB/experimentalSprayData/', campaignID]);
            end
            
            disp(['    Saving to: ', saveLoc, '/Experimental/MATLAB/experimentalSprayData/', campaignID, '/', caseID, '.mat']);
            save([saveLoc, '/Experimental/MATLAB/experimentalSprayData/', campaignID, '/', caseID, '.mat'], ...
                 'campaignID', 'caseID', 'planeID', 'expSprayData', 'sampleFreq', '-v7.3', '-noCompression');
            disp('        Success');
            
            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end
        
    end
    
end