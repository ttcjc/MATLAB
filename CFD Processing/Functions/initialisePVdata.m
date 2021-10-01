%% ParaView Data Initialisation v1.0
% ----
% Collates Basic Case Data for Further Processing of ParaView Exports
% ----
% Usage: [caseFolder, xDims, yDims, zDims, geometry] = initialisePVdata(field)
%        'field' -> Desired Field Stored as String


%% Changelog

% v1.0 - Initial Commit


%% Supported Case Types

% Lag_Test (No Geometry)
% Test_Block
% Windsor (Balance)
% Windsor (Upstream)


%% Main Function

function [caseFolder, data, xDims, yDims, zDims, geometry] = initialisePVdata(field)

    disp('CASE SELECTION');
    disp('--------------');
    disp(' ');

    caseFolder = uigetdir('~/Documents/Engineering/PhD/Data/Numerical/OpenFOAM', ...
                          'Select Case');

    disp(['Case: ', caseFolder]);
    disp(' ');
    
    % Confirm Support
    if ~contains(caseFolder, ["Lag_Test", "Test_Block", "Windsor"])
        error('Invalid Case Directory (Unsupported Case Type)');
    end
    
    % Confirm Data Availability
    dataFiles = dir([caseFolder, '/*.csv']);

    if isempty(dataFiles)
        error('Invalid Case Directory (No Data Files Available)');
    end

    j = 1;
    for i = 1:height(dataFiles)

        if contains(dataFiles(i,1).name, [field, '_'])
            data.files{j,1} = dataFiles(i,1).name;
            j = j + 1;
        end

    end

    if isempty(data.files)
        error('Invalid Case Directory (No Valid Data Available)');
    end

    % Set Dimensions of Interest
    if contains(caseFolder, 'Lag_Test')
        xDims = [-14; 14];
        yDims = [-0.96; 0.96];
        zDims = [0; 1.32];
    elseif contains(caseFolder, 'Test_Block')
        xDims = [-0.56075; 0.48325];
        yDims = [-0.1945; 0.1945];
        zDims = [0.05; 0.339];
    elseif contains(caseFolder, 'Windsor') && contains(caseFolder, 'Balance')
        xDims = [-0.56075; 0.48325];
        yDims = [-0.1945; 0.1945];
        zDims = [0.05; 0.339];
    elseif contains(caseFolder, 'Windsor') && contains(caseFolder, 'Upstream')
        xDims = [-1.88575; -0.84175];
        yDims = [-0.1945; 0.1945];
        zDims = [0.05; 0.339];
    end
    
    disp(' ');

    disp('GEOMETRY SELECTION');
    disp('------------------');
    disp(' ');

    if contains(caseFolder, 'Lag_Test')
        disp('No Test Geometry Required');
        geometry = 0;
    else
        [file, path] = uigetfile('~/CAD/CFD Geometries/*.stl', ...
                                 'Select Subject Geometry', 'multiSelect', 'on');

        if isa(file, 'cell')
            disp('Geometries: ');

            for i = 1:width(file)
                part = file{1,i}(1:end-4);
                geometry.(part) = stlreader([path,file{1,i}]);
                disp(['    ', part]);
                
                if contains(part, 'Windsor')
                    geometry.(part).vertices = geometry.(part).vertices / 1.044;
                end
                
            end

        else
            part = file(1:end-4);
            geometry.(part) = stlreader([path, file]);
            disp(['Geometry: ', part]);
            
            if contains(part, 'Windsor')
                geometry.(part).vertices = geometry.(part).vertices / 1.044;
            end
            
        end
    end

end