%% Volume Field Plotter v2.6
% ----
% Plots Previously Processed Volume Fields
% ----
% Usage: [fig, surfaceNo] = plotVolumeField(xLimsData, yLimsData, zLimsData, spatialRes, ...
%                                           xOrig, yOrig, zOrig, POD, fieldData, nSurfaces, surfaceNo, ...
%                                           fig, figName, geometry, isoValue, cMap, figTitle, viewAngle, ...
%                                           multiView, xLimsPlot, yLimsPlot, zLimsPlot, figSave);
%
%        '*LimsData'   -> Contour Plot Limits [Dimensions of '*Orig']
%        'spatialRes'  -> Target Grid Spacing [Dimensions of '*Orig']
%        '*Orig'       -> Origial 3D Arrays of Cartesian Positions
%        'POD'         -> POD Mode Presentation [True/False]
%        'fieldData'   -> 3D Array of Field Data @ '*Orig' Points
%        'nSurfaces'   -> Number of Isosurfaces in a Multi-Surface Figure
%        'surfaceNo'   -> Current Surface Number
%        'fig'         -> Figure Number
%        'figName'     -> Figure Name
%        'geometry'    -> STL(s) to Include in Plot
%        'isoValue'    -> Field Value Used for Isosurface
%        'cMap'        -> Colour Map
%        'figTitle'    -> Figure Title
%        'viewAngle'   -> Default Viewing Angle
%        'multiView'   -> Plot Additional Views Normal to Each Co-Ordinate Axis [True/False]
%        '*LimsPlot'   -> 3D Axes Limits [Dimensions of '*Orig']
%        'figSave'     -> Save .fig File [True/False]


%% Changelog

% v1.0 - Origial Commit
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats
% v2.1 - Cleaned-up POD Functionality
% v2.2 - Rename and Minor Formatting Updates
% v2.3 - Added Spatial Resolution as an Input Variable
% v2.4 - Added Support for Saving Multiple View Angles
% v2.5 - Update To Ensure Consistent Figure Sizes
% v2.6 - Added Support for Multiple Surfaces


%% Main Function

function [fig, surfaceNo] = plotVolumeField(xLimsData, yLimsData, zLimsData, spatialRes, ...
                                            xOrig, yOrig, zOrig, POD, fieldData, nSurfaces, surfaceNo, ...
                                            fig, figName, geometry, isoValue, cMap, figTitle, viewAngle, ...
                                            multiView, xLimsPlot, yLimsPlot, zLimsPlot, figSave)
    
    % Generate Refined Grid
    cellSizeX = (xLimsData(2) - xLimsData(1)) / round(((xLimsData(2) - xLimsData(1)) / spatialRes));
    cellSizeY = (yLimsData(2) - yLimsData(1)) / round(((yLimsData(2) - yLimsData(1)) / spatialRes));
    cellSizeZ = (zLimsData(2) - zLimsData(1)) / round(((zLimsData(2) - zLimsData(1)) / spatialRes));
    
    [x, y, z] = ndgrid(xLimsData(1):cellSizeX:xLimsData(2), ...
                       yLimsData(1):cellSizeY:yLimsData(2), ...
                       zLimsData(1):cellSizeZ:zLimsData(2));
    
    % Convert From 'ndgrid' to 'meshgrid' Format
    xOrig = permute(xOrig, [2,1,3]);
    yOrig = permute(yOrig, [2,1,3]);
    zOrig = permute(zOrig, [2,1,3]);
    x = permute(x, [2,1,3]);
    y = permute(y, [2,1,3]);
    z = permute(z, [2,1,3]);

    % Smooth Data
    if POD
        
        for i = 1:height(fieldData)
            fieldData{i} = permute(fieldData{i}, [2,1,3]);
            fieldData{i} = interp3(xOrig, yOrig, zOrig, fieldData{i}, x, y, z);
            fieldData{i} = smooth3(fieldData{i}, 'box', 3);
        end
        
    else
        fieldData = permute(fieldData, [2,1,3]);
        fieldData = interp3(xOrig, yOrig, zOrig, fieldData, x, y, z);
        fieldData = smooth3(fieldData, 'box', 3);
    end
    
    % Origialise Figure
    if surfaceNo == 1
        fig = fig + 1;
        set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                         'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
        pause(0.5);
        hold on;
        set(gca, 'positionConstraint', 'outerPosition', 'dataAspectRatio', [1, 1, 1], ...
                 'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');
        lighting gouraud;
        
        % Plot Geometry
        if ~isempty(geometry)
            parts = fieldnames(geometry);
            
            for i = 1:height(parts)
                patch('faces', geometry.(parts{i}).faces, ...
                      'vertices', geometry.(parts{i}).vertices, ...
                      'faceColor', [0.5, 0.5, 0.5], ...
                      'edgeColor', [0.5, 0.5, 0.5], ...
                      'lineStyle', 'none');
            end
            clear i;
            
        end
        
    end

    % Plot Iso-Surface
    if POD
        iso = isosurface(x, y, z, fieldData{1}, isoValue);
        p = patch(iso, 'faceColor', cMap(1,:), 'edgeColor', 'none');
        isonormals(fieldData{1}, p);
        
        iso = isosurface(x, y, z, fieldData{2}, isoValue);
        p = patch(iso, 'faceColor', cMap(end,:), 'edgeColor', 'none');
        isonormals(fieldData{2}, p);
    else
        iso = isosurface(x, y, z, fieldData, isoValue);
        p = patch(iso, 'faceColor', cMap, 'edgeColor', 'none');
        isonormals(fieldData, p);
    end
    
    % Format Figure
    if surfaceNo == nSurfaces
        title('{-----}', 'interpreter', 'latex');
        subtitle(figTitle);
        lightangle(0, 45);
        axis on;
        box on;
        grid off;
        view(viewAngle);
        xlim([xLimsPlot(1), xLimsPlot(2)]);
        ylim([yLimsPlot(1), yLimsPlot(2)]);
        zlim([zLimsPlot(1), zLimsPlot(2)]);
        tickData = [];
        xticks(tickData);
        tickData = [];
        yticks(tickData);
        tickData = [];
        zticks(tickData);
        hold off;

        tightInset = get(gca, 'TightInset');
        set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                                   (tightInset(2) + 0.00625), ...
                                   (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                                   (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
        pause(0.5);
        hold off;

        % Save Figure
        if multiView
            multiAngle = [viewAngle; 0, 0; 0, 90; 90, 0];
            
            for i = 1:height(multiAngle)
                view(multiAngle(i,:));
                
                print(gcf, [userpath, '/Output/Figures/', figName, '_', num2str(i), '.png'], '-dpng', '-r300');
                
                pause(0.5);
                
                if figSave
                    savefig(gcf, [userpath, '/Output/Figures/', figName, '_', num2str(i), '.fig']);
                end
                
            end
            
        else
            print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');
            
            if figSave
                savefig(gcf, [userpath, '/Output/Figures/', figName, '.fig']);
            end
            
        end
        
    else
        surfaceNo = surfaceNo + 1;
    end
    
end
