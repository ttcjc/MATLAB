%% Geometry Selection Tool v1.0
% ----
% Load One or More ASCII STL Geometry Files for Further Processing
% ----
% Usage: [geometry, xDims, yDims, zDims, precision, normalise] = selectGeometry(normalise);
%        'normalise' -> Normalise Dimensions [True/False]


%% Changelog

% v1.0 - Initial Commit


%% Supported Geometries

% Test Block
% Windsor Model


%% Main Function

function [geometry, xDims, yDims, zDims, spacePrecision, normalise] = selectGeometry(normalise)

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

    disp(' ');

    disp('Defining Geometry Boundaries...');
    
    % Define Geometric Boundaries
    xDims = nan(2,1);
    yDims = xDims;
    zDims = xDims;
    
    parts = fieldnames(geometry);
    
    for i = 1:height(parts)
        geometry.(parts{i}).boundaries.YZ = cell(height(parts),1);
        geometry.(parts{i}).boundaries.XZ = geometry.(parts{i}).boundaries.YZ;
        geometry.(parts{i}).boundaries.XY = geometry.(parts{i}).boundaries.YZ;

        geoPoints = geometry.(parts{i}).vertices;

        % 3D Boundaries
        index = boundary(geoPoints(:,1), geoPoints(:,2), geoPoints(:,3), 1);
        geoPoints = geoPoints(index,:);
        
        xDims(1) = min(xDims(1), min(geoPoints(:,1)));
        xDims(2) = max(xDims(2), max(geoPoints(:,1)));
        yDims(1) = min(yDims(1), min(geoPoints(:,2)));
        yDims(2) = max(yDims(2), max(geoPoints(:,2)));
        zDims(1) = min(zDims(1), min(geoPoints(:,3)));
        zDims(2) = max(zDims(2), max(geoPoints(:,3)));
        
        % YZ Boundaries
        index = boundary(geoPoints(:,2), geoPoints(:,3), 0);
        geometry.(parts{i}).boundaries.YZ = nan(height(index),3);
        geometry.(parts{i}).boundaries.YZ(:,[2,3]) = geoPoints(index,[2,3]);
        
        % XZ Boundaries
        index = boundary(geoPoints(:,1), geoPoints(:,3), 0);
        geometry.(parts{i}).boundaries.XZ = nan(height(index),3);
        geometry.(parts{i}).boundaries.XZ(:,[1,3]) = geoPoints(index,[1,3]);
        
        % XY Boundaries
        index = boundary(geoPoints(:,1), geoPoints(:,2), 0);
        geometry.(parts{i}).boundaries.XY = nan(height(index),3);
        geometry.(parts{i}).boundaries.XY(:,[1,2]) = geoPoints(index,[1,2]);
    end
    
    % Define Rounding Precision and Normalise Dimensions
    xPre = max(width(extractAfter(num2str(xDims(1), 7), '.')), width(extractAfter(num2str(xDims(2), 7), '.')));
    yPre = max(width(extractAfter(num2str(yDims(1), 7), '.')), width(extractAfter(num2str(yDims(2), 7), '.')));
    zPre = max(width(extractAfter(num2str(zDims(1), 7), '.')), width(extractAfter(num2str(zDims(2), 7), '.')));
    
    spacePrecision = max([xPre, yPre, zPre]);

    if normalise
        disp(' ');
        
        disp('Normalising Dimensions...');
        
        if contains(path, ["Run_Test", "Windsor"])
            xDims = round((xDims / 1.044), spacePrecision);
            yDims = round((yDims / 1.044), spacePrecision);
            zDims = round((zDims / 1.044), spacePrecision);

            for i = 1:height(parts)
                geometry.(parts{i}).vertices = round((geometry.(parts{i}).vertices / 1.044), spacePrecision);
                geometry.(parts{i}).boundaries.YZ(:,[2,3]) = round((geometry.(parts{i}).boundaries.YZ(:,[2,3]) / 1.044), spacePrecision);
                geometry.(parts{i}).boundaries.XZ(:,[1,3]) = round((geometry.(parts{i}).boundaries.XZ(:,[1,3]) / 1.044), spacePrecision);
                geometry.(parts{i}).boundaries.XY(:,[1,2]) = round((geometry.(parts{i}).boundaries.XY(:,[1,2]) / 1.044), spacePrecision);
            end

        else
            normalise = false;
            disp('    WARNING: Dimension Normalisation for This Case Type Not Supported');
        end
    
    end

end