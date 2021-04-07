%% Case Initialisation v1.1
% ----
% Collates Basic Case Data for Further Processing
% ----
% Usage: [caseFolder, xDims, yDims, zDims, timeDirs, deltaT, geometry] = initialiseCase

%% Changelog

% v1.0 - Initial Commit
% v1.1 - Added Support for Balance and Upstream Windsor Case Variants


%% Supported Case Types

% Lag_Test (No Geometry)
% Test_Block
% Windsor (Balance)
% Windsor (Upstream)


%% Main Function

function [caseFolder, xDims, yDims, zDims, timeDirs, deltaT, geometry] = initialiseCase

    disp('CASE SELECTION');
    disp('--------------');
    disp(' ');

    caseFolder = uigetdir('~/OpenFOAM', 'Select Case');

    disp(['Case: ', caseFolder]);
    disp(' ');
    disp(' ');

    % Confirm Case Validity and Identify Time Directories
    [timeDirs, deltaT] = timeDirectories(caseFolder);
    
    % Confirm Support
    if ~contains(caseFolder, ["Lag_Test", "Test_Block", "Windsor"])
        error('Invalid Case Directory (Unsupported Case Type)');
    end

    % Set Dimensions of Interest
    if contains(caseFolder, 'Lag_Test')
        xDims = [-14; 14];
        yDims = [-0.96; 0.96];
        zDims = [0; 1.32];
    elseif contains(caseFolder, 'Test_Block')
        xDims = [-0.56075; 0.48325];
        yDims = [-(0.389 / 2); (0.389 / 2)];
        zDims = [0.05; 0.339];
    elseif contains(caseFolder, 'Windsor') && contains(caseFolder, 'Balance')
        xDims = [-0.56075; 0.48325];
        yDims = [-(0.389 / 2); (0.389 / 2)];
        zDims = [0.05; 0.339];
    elseif contains(caseFolder, 'Windsor') && contains(caseFolder, 'Upstream')
        xDims = [-1.88575; -0.84175];
        yDims = [-(0.389 / 2); (0.389 / 2)];
        zDims = [0.05; 0.339];
    end
    
    disp(' ');
    disp(' ');

    disp('GEOMETRY SELECTION');
    disp('------------------');
    disp(' ');

    if contains(caseFolder, 'Lag_Test')
        disp('No Test Geometry Required');
        geometry = 0;
    else
        [file, path] = uigetfile('~/CAD/CFD Geometries/*.stl', 'Select Subject Geometry', 'multiSelect', 'on');

        if isa(file, 'cell')
            disp('Geometries: ');

            for i = 1:size(file,2)
                part = file{1,i}(1:end-4);
                geometry.(part) = stlreader([path,file{1,i}]);
                disp(['    ', part]);
            end

        else
            part = file(1:end-4);
            geometry.(part) = stlreader([path, file]);
            disp(['Geometry: ', part]);
        end
    end

end