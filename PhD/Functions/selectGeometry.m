%% Geometry Selection Tool v1.0
% ----
% Load One or More ASCII STL Geometry Files for Further Processing
% ----
% Usage: [geometry, xDims, yDims, zDims, precision] = selectGeometry(normalise);
%        'normalise' -> Normalise Dimensions [True/False]


%% Changelog

% v1.0 - Initial Commit


%% Supported Geometries

% Test Block
% Windsor Model


%% Main Function

function [geometry, xDims, yDims, zDims, precision] = selectGeometry(normalise)

    disp('Geometry Selection');
    disp('-------------------');
    
    disp(' ');
    
    % Select STL File(s)
    [file, path] = uigetfile('~/CAD/CFD Geometries/*.stl', ...
                             'Select Subject Geometry', 'multiSelect', 'on');
    
    if isa(file, 'cell')
        disp('Loading Geometries:');
        
        for i = 1:width(file)
            parts = file{i}(1:end-4); % Ignore '.stl'
            geometry.(parts) = stlreader([path,file{i}]);
            disp(['    ', parts]);
        end
        
    else
        parts = file(1:end-4); % Ignore '.stl'
        geometry.(parts) = stlreader([path, file]);
        disp('Loading Geometry:');
        disp(['    ', parts]);
    end

    disp (' ');

    disp('Defining Bounding Box...');

    % Define Bounding Box
    if contains(path, 'Run_Test')
        xDims = [-1.88575; -0.84175];
        yDims = [-0.1945; 0.1945];
        zDims = [0.05; 0.339];
    elseif contains(path, 'Windsor') && contains(path, 'Balance')
        xDims = [-0.56075; 0.48325];
        yDims = [-0.1945; 0.1945];
        zDims = [0.05; 0.339];
    elseif contains(path, 'Windsor') && contains(path, 'Upstream')
        xDims = [-1.88575; -0.84175];
        yDims = [-0.1945; 0.1945];
        zDims = [0.05; 0.339];
    else
        error('Invalid Geometry (Unsupported Case Type)');

    end
    
    % Define Rounding Precision and Normalise Dimensions
    xPre = max(width(extractAfter(num2str(xDims(1), 8), '.')), width(extractAfter(num2str(xDims(2), 8), '.')));
    yPre = max(width(extractAfter(num2str(yDims(1), 8), '.')), width(extractAfter(num2str(yDims(2), 8), '.')));
    zPre = max(width(extractAfter(num2str(zDims(1), 8), '.')), width(extractAfter(num2str(zDims(2), 8), '.')));
    
    precision = max([xPre, yPre, zPre]);

    if normalise
        disp (' ');
        
        disp('Normalising Dimensions...');
        
        if contains(path, ["Run_Test", "Windsor"])
            xDims = round(xDims / 1.044, precision);
            yDims = round(yDims / 1.044, precision);
            zDims = round(zDims / 1.044, precision);

            parts = fieldnames(geometry);

            for i = 1:height(parts)
                geometry.(parts{i}).vertices = round(geometry.(parts{i}).vertices / 1.044, precision);
            end

        else
            disp('    WARNING: Dimension Normalisation for This Case Type Not Supported');
        end
    
    end

end