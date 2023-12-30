%% Geometry Selection Tool v1.4
% ----
% Load One or More ASCII STL Geometry Files for Further Processing
% ----
% Usage: [geometry, xDims, yDims, zDims, precision, normLength] = selectGeometry(geoLoc);
% 
%        'geoLoc' -> Start of File Path, Stored as a String


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Formatting Updates
% v1.2 - Added Support for Full-Scale Windsor Model Geometries
% v1.3 - Fixed a Bug Preventing the Normalisation of the Quarter-Scale Windsor Model
% v1.4 - Removed Geometry Normalisation (Should Be Done After Processing)


%% Supported Geometries

% Test Block
% Windsor Model (Quarter-Scale & Full-Scale)


%% Main Function

function [geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry(geoLoc)

    disp('Geometry Selection');
    disp('-------------------');
    
    disp(' ');
    
    % Select STL File(s)
    [file, path] = uigetfile([geoLoc, '/*.stl'], ...
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

    disp(' ');

    disp('Defining Geometry Boundaries...');
    
    % Define Geometric Boundaries
    xDims = nan(2,1);
    yDims = xDims;
    zDims = xDims;
    
    parts = fieldnames(geometry);
    
    for i = 1:height(parts)
        geoPoints = geometry.(parts{i}).vertices;
        index = boundary(geoPoints(:,1), geoPoints(:,2), geoPoints(:,3), 1);
        geoPoints = geoPoints(index,:);
        
        xDims(1) = min(xDims(1), min(geoPoints(:,1)));
        xDims(2) = max(xDims(2), max(geoPoints(:,1)));
        yDims(1) = min(yDims(1), min(geoPoints(:,2)));
        yDims(2) = max(yDims(2), max(geoPoints(:,2)));
        zDims(1) = min(zDims(1), min(geoPoints(:,3)));
        zDims(2) = max(zDims(2), max(geoPoints(:,3)));
    end
    
    disp('    Success');
    
    % Define Rounding Precision and Normalise Dimensions
    xPre = max(width(extractAfter(num2str(xDims(1), 7), '.')), width(extractAfter(num2str(xDims(2), 7), '.')));
    yPre = max(width(extractAfter(num2str(yDims(1), 7), '.')), width(extractAfter(num2str(yDims(2), 7), '.')));
    zPre = max(width(extractAfter(num2str(zDims(1), 7), '.')), width(extractAfter(num2str(zDims(2), 7), '.')));
    
    spacePrecision = max([xPre, yPre, zPre]);
    
    % Specify Value Used for Dimensional Normalisation
    normLength = round((diff(xDims)), spacePrecision);
    
end