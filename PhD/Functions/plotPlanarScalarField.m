%% Planar Scalar Field Plotter v2.1
% ----
% Plots Previously Processed Planar Scalar Fields
% ----
% Usage: [fig, planeNo] = plotPlanarScalarField(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
%                                               mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, contourlines, ...
%                                               xDims, yDims, zDims, refPoint, figTitle, figSubtitle, cLims, ...
%                                               xLimsPlot, yLimsPlot, zLimsPlot, normalise, nPlanes, planeNo);
%        'orientation'  -> Plane Orientation ['YZ', 'XZ', 'XY']
%        '*LimsData'    -> Contour Plot Limits
%        'positionData' -> Cartesian Positions of Data Points
%        'scalarData'   -> Field Value @ 'positionData' Points
%        'mapPerim'     -> Map Perimeter Used When Plotting a Plane of Arbitrary Shape
%        'nPlanes'      -> Number of Planes in a Multi-Plane Figure
%        'planeNo'      -> Current Plane Number
%        'fig'          -> Figure Number
%        'figName'      -> Figure Name
%        'cMap'         -> Colour Map
%        'geometry'     -> STL(s) to Include in Plot
%        'contourlines' -> Contour Lines at the Corresponding Values
%        '*Dims'        -> Simple Bounding Box of Geometry
%        'refPoint'     -> Reference Point (CoM, CoP etc.)
%        'figTitle'     -> Leave Blank ('-') for Formatting Purposes
%        'figSubtitle'  -> Figure Title
%        'cLims'        -> Colour Map Limits
%        '*LimsPlot'    -> 3D Axes Limits
%        'normalise'    -> Normalise Dimensions [True/False]


%% Changelog

% v1.0 - Initial Commit (Expanded Functionality of 'contaminantPlots')
% v2.0 - Shifted to Using 'griddedInterpolant' and Added Support for Multi-Plane Figures)
% v2.1 - Rename and Minor Formatting Updates


%% Main Function

function [fig, planeNo] = plotPlanarScalarField(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                                mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, contourlines, ...
                                                xDims, yDims, zDims, refPoint, figTitle, figSubtitle, cLims, ...
                                                xLimsPlot, yLimsPlot, zLimsPlot, normalise)
    
    cellSize = 1e-3; % [m or l]
    
    % Format Data
    switch orientation
        
        case 'YZ'
            % Reshape Data for Improved Interpolation Performance
            gridShape = [height(unique(positionData(:,2))), ...
                         height(unique(positionData(:,3)))];

            y = reshape(positionData(:,2), gridShape);
            z = reshape(positionData(:,3), gridShape);
            
            scalar = reshape(scalarData, gridShape);
            
            % Perform Interpolation
            interp = griddedInterpolant(y, z, scalar, 'linear', 'none');
            
            % Generate 3D Surface
            cellSizeX = cellSize;
            cellSizeY = (yLimsData(2) - yLimsData(1)) / round(((yLimsData(2) - yLimsData(1)) / cellSize));
            cellSizeZ = (zLimsData(2) - zLimsData(1)) / round(((zLimsData(2) - zLimsData(1)) / cellSize));
            
            [x, y, z] = ndgrid((xLimsData - cellSizeX):cellSizeX:(xLimsData + cellSizeX), ...
                               yLimsData(1):cellSizeY:yLimsData(2), ...
                               zLimsData(1):cellSizeZ:zLimsData(2));
            
            % Map Data on to 3D Surface
            scalar = zeros(size(x));
            scalar(2,:,:) = interp(y(2,:,:), z(2,:,:));
            
            % Remove Data Outside Desired Map Perimeter
            if ~isempty(mapPerim)
                [indexIn, indexOn] = inpolygon(y, z, mapPerim(:,2), mapPerim(:,3));
                indexMap = double(or(indexIn, indexOn));
                indexMap(indexMap == 0) = nan;
                
                scalar = scalar .* indexMap;
            end
            
        case 'XZ'
            % Reshape Data for Improved Interpolation Performance
            gridShape = [height(unique(positionData(:,1))), ...
                         height(unique(positionData(:,3)))];

            x = reshape(positionData(:,1), gridShape);
            z = reshape(positionData(:,3), gridShape);
            
            scalar = reshape(scalarData, gridShape);
            
            % Perform Interpolation
            interp = griddedInterpolant(x, z, scalar, 'linear', 'none');
            
            % Generate 3D Surface
            cellSizeX = (xLimsData(2) - xLimsData(1)) / round(((xLimsData(2) - xLimsData(1)) / cellSize));
            cellSizeY = cellSize;
            cellSizeZ = (zLimsData(2) - zLimsData(1)) / round(((zLimsData(2) - zLimsData(1)) / cellSize));

            [x, y, z] = ndgrid(xLimsData(1):cellSizeX:xLimsData(2), ...
                               (yLimsData - cellSizeY):cellSizeY:(yLimsData + cellSizeY), ...
                               zLimsData(1):cellSizeZ:zLimsData(2));
    
            % Map Scalar Data on to 3D Surface
            scalar = zeros(size(x));
            scalar(:,2,:) = interp(x(:,2,:), z(:,2,:));
            
            % Remove Data Outside Desired Map Perimeter
            if ~isempty(mapPerim)
                [indexIn, indexOn] = inpolygon(x, z, mapPerim(:,1), mapPerim(:,3));
                indexMap = double(or(indexIn, indexOn));
                indexMap(indexMap == 0) = nan;
                
                scalar = scalar .* indexMap;
            end
            
        case 'XY'
            % Reshape Data for Improved Interpolation Performance
            gridShape = [height(unique(positionData(:,1))), ...
                         height(unique(positionData(:,2)))];

            x = reshape(positionData(:,1), gridShape);
            y = reshape(positionData(:,2), gridShape);
            
            scalar = reshape(scalarData, gridShape);
            
            % Perform Interpolation
            interp = griddedInterpolant(x, y, scalar, 'linear', 'none');
            
            % Generate 3D Surface
            cellSizeX = (xLimsData(2) - xLimsData(1)) / round(((xLimsData(2) - xLimsData(1)) / cellSize));
            cellSizeY = (yLimsData(2) - yLimsData(1)) / round(((yLimsData(2) - yLimsData(1)) / cellSize));
            cellSizeZ = cellSize;
            
            [x, y, z] = ndgrid(xLimsData(1):cellSizeX:xLimsData(2), ...
                               yLimsData(1):cellSizeY:yLimsData(2), ...
                               (zLimsData - cellSizeZ):cellSizeZ:(zLimsData + cellSizeZ));

            % Map Scalar Data on to 3D Surface
            scalar = zeros(size(x));
            scalar(:,:,2) = interp(x(:,:,2), y(:,:,2));
            
            % Remove Data Outside Desired Map Perimeter
            if ~isempty(mapPerim)
                [indexIn, indexOn] = inpolygon(x, y, mapPerim(:,1), mapPerim(:,2));
                indexMap = double(or(indexIn, indexOn));
                indexMap(indexMap == 0) = nan;
                
                scalar = scalar .* indexMap;
            end
            
    end

    % Present Data
    switch orientation
        
        case 'YZ'
            % Initialise Figure
            if planeNo == 1
                fig = fig + 1;
%                 set(figure(fig), 'color', [1, 1, 1], 'outerPosition', [25, 25, 850, 850], 'name', figName);
                set(figure(fig), 'name', figName, 'color', [1, 1, 1], 'paperPositionMode', 'manual', 'paperUnits', 'inches', ...
                                 'paperSize', [3.45, 3.45], 'paperPosition', [0.05, 0.05, 3.35, 3.35]);
                set(gca, 'dataAspectRatio', [1, 1, 1], 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
                         'fontSize', 16, 'layer', 'top');
                lighting gouraud;
                colormap(cMap);
                hold on;
                
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
                    
                end

            end
            
            % Plot Contour
            surf(squeeze(x(2,:,:)), squeeze(y(2,:,:)), squeeze(z(2,:,:)), squeeze(scalar(2,:,:)), ...
                 'lineStyle', 'none', 'faceLighting', 'none');

            if nPlanes == 1 && ~isempty(geometry) && (xLimsData < xDims(1) || xLimsData > xDims(2))
                
                for i = 1:height(parts)
                    geometry.(parts{i}).boundaries.YZ(:,1) = xLimsData;
                    
                    plot3(geometry.(parts{i}).boundaries.YZ(:,1), ...
                          geometry.(parts{i}).boundaries.YZ(:,2), ...
                          geometry.(parts{i}).boundaries.YZ(:,3), ...
                          'color', 'w', 'lineStyle', '-', 'lineWidth', 1.5);
                end
                
            end
            
            % Plot Contour Lines
            if ~isempty(contourlines)
                % Convert From 'ndgrid' to 'meshgrid' Format
                x = permute(x, [2,1,3]);
                y = permute(y, [2,1,3]);
                z = permute(z, [2,1,3]);
                scalar = permute(scalar, [2,1,3]);
                
                contours = contourslice(x, y, z, scalar, xLimsData, [], [], contourlines);
                set(contours, 'edgeColor', 'w', 'lineStyle', '-', 'lineWidth', 1.5);
            end
            
            % Plot Reference Point
            if ~isempty(refPoint)
                scatter3(refPoint(1), refPoint(2), refPoint(3), 200, 'w', 'x', 'lineWidth', 2);
            end
            
            % Format Figure
            if planeNo == nPlanes

                if nPlanes == 1
                    title(figTitle, 'color', ([254, 254, 254] / 255));
                    subtitle(figSubtitle);
                    lightangle(90, 45);
                    axis on;
                    box on;
                    grid off;
                    caxis(cLims);
                    view([90, 0]);
                    xlim([xLimsPlot(1), xLimsPlot(2)]);
                    ylim([yLimsPlot(1), yLimsPlot(2)]);
                    zlim([zLimsPlot(1), zLimsPlot(2)]);
                    tickData = [];
                    xticks(tickData);
                    tickData = round((yLimsPlot(1):((yLimsPlot(2) - yLimsPlot(1)) / 5):yLimsPlot(2)), 2);
                    yticks(tickData(2:(end - 1)));
                    tickData = round((zLimsPlot(1):((zLimsPlot(2) - zLimsPlot(1)) / 5):zLimsPlot(2)), 2);
                    zticks(tickData(2:(end - 1)));
                    xtickformat('%+.2f');
                    ytickformat('%+.2f');
                    ztickformat('%+.2f');
                    
                    if normalise
                        ylabel('{\bf{y}}_{{\it{l}}}', 'fontName', 'LM Roman 12');
                        zlabel('{\bf{z}}_{{\it{l}}}', 'fontName', 'LM Roman 12');
                    else
                        ylabel('{\bf{y}}_{{\it{m}}}', 'fontName', 'LM Roman 12');
                        zlabel('{\bf{z}}_{{\it{m}}}', 'fontName', 'LM Roman 12');
                    end
                
                else
                    title(figTitle, 'color', ([254, 254, 254] / 255));
                    subtitle(figSubtitle);
                    lightangle(0, 45);
                    axis on;
                    box on;
                    grid off;
                    view([30, 30]);
                    xlim([xLimsPlot(1), xLimsPlot(2)]);
                    ylim([yLimsPlot(1), yLimsPlot(2)]);
                    zlim([zLimsPlot(1), zLimsPlot(2)]);
                    tickData = [];
                    xticks(tickData);
                    tickData = [];
                    yticks(tickData);
                    tickData = [];
                    zticks(tickData);
                end
                
%                 set(gca, 'outerPosition', [0.05, 0.05, 0.9, 0.9]);
                hold off;
                
                % Save Figure
                pause(2);
                
                exportgraphics(gcf, [userpath, '/Output/Figures/', figName, '.png'], 'resolution', 600);
%                 print(gcf, [userpath, '/Output/Figures/', figName, '.pdf'], '-image', '-dpdf', '-r600')
                savefig(gcf, [userpath, '/Output/Figures/', figName, '.fig']);
            else
                planeNo = planeNo + 1;
            end
            
        case 'XZ'
            % Initialise Figure
            if planeNo == 1
                fig = fig + 1;
%                 set(figure(fig), 'color', [1, 1, 1], 'outerPosition', [25, 25, 850, 850], 'name', figName);
                set(figure(fig), 'name', figName, 'color', [1, 1, 1], 'paperPositionMode', 'manual', 'paperUnits', 'inches', ...
                                 'paperSize', [3.45, 3.45], 'paperPosition', [0.05, 0.05, 3.35, 3.35]);
                set(gca, 'dataAspectRatio', [1, 1, 1], 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
                         'fontSize', 16, 'layer', 'top');
                lighting gouraud;
                colormap(cMap);
                hold on;
                
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
                    
                end

            end
            
            % Plot Contour
            surf(squeeze(x(:,2,:)), squeeze(y(:,2,:)), squeeze(z(:,2,:)), squeeze(scalar(:,2,:)), ...
                 'lineStyle', 'none', 'faceLighting', 'none');

            if nPlanes == 1 && ~isempty(geometry) && yLimsData < yDims(1)
                
                for i = 1:height(parts)
                    geometry.(parts{i}).boundaries.XZ(:,2) = yLimsData;
                    
                    plot3(geometry.(parts{i}).boundaries.XZ(:,1), ...
                          geometry.(parts{i}).boundaries.XZ(:,2), ...
                          geometry.(parts{i}).boundaries.XZ(:,3), ...
                          'color', 'w', 'lineStyle', '-', 'lineWidth', 1.5);
                end
                
            end
            
            % Plot Contour Lines
            if ~isempty(contourlines)
                % Convert From 'ndgrid' to 'meshgrid' Format
                x = permute(x, [2,1,3]);
                y = permute(y, [2,1,3]);
                z = permute(z, [2,1,3]);
                scalar = permute(scalar, [2,1,3]);
                
                contours = contourslice(x, y, z, scalar, [], yLimsData, [], contourlines);
                set(contours, 'edgeColor', 'w', 'lineStyle', '-', 'lineWidth', 1.5);
            end
            
            % Plot Reference Point
            if ~isempty(refPoint)
                scatter3(refPoint(1), refPoint(2), refPoint(3), 200, 'w', 'x', 'lineWidth', 2);
            end
            
            % Format Figure
            if planeNo == nPlanes

                if nPlanes == 1
                    title(figTitle, 'color', ([254, 254, 254] / 255));
                    subtitle(figSubtitle);
                    lightangle(0, 45);
                    axis on;
                    box on;
                    grid off
                    caxis(cLims);
                    view([0, 0]);
                    xlim([xLimsPlot(1), xLimsPlot(2)]);
                    ylim([yLimsPlot(1), yLimsPlot(2)]);
                    zlim([zLimsPlot(1), zLimsPlot(2)]);
                    tickData = round((xLimsPlot(1):((xLimsPlot(2) - xLimsPlot(1)) / 5):xLimsPlot(2)), 2);
                    xticks(tickData(2:(end - 1)));
                    tickData = [];
                    yticks(tickData);
                    tickData = round((zLimsPlot(1):((zLimsPlot(2) - zLimsPlot(1)) / 5):zLimsPlot(2)), 2);
                    zticks(tickData(2:(end - 1)));
                    xtickformat('%+.2f');
                    ytickformat('%+.2f');
                    ztickformat('%+.2f');
                    
                    if normalise
                        xlabel('{\bf{x}}_{{\it{l}}}', 'fontName', 'LM Roman 12');
                        zlabel('{\bf{z}}_{{\it{l}}}', 'fontName', 'LM Roman 12');
                    else
                        xlabel('{\bf{x}}_{{\it{m}}}', 'fontName', 'LM Roman 12');
                        zlabel('{\bf{z}}_{{\it{m}}}', 'fontName', 'LM Roman 12');
                    end

                else
                    title(figTitle, 'color', ([254, 254, 254] / 255));
                    subtitle(figSubtitle);
                    lightangle(0, 45);
                    axis on;
                    box on;
                    grid off;
                    view([30, 30]);
                    xlim([xLimsPlot(1), xLimsPlot(2)]);
                    ylim([yLimsPlot(1), yLimsPlot(2)]);
                    zlim([zLimsPlot(1), zLimsPlot(2)]);
                    tickData = [];
                    xticks(tickData);
                    tickData = [];
                    yticks(tickData);
                    tickData = [];
                    zticks(tickData);
                end
                
%                 set(gca, 'outerPosition', [0.05, 0.05, 0.9, 0.9]);
                hold off;
                
                % Save Figure
                pause(2);
                
                exportgraphics(gcf, [userpath, '/Output/Figures/', figName, '.png'], 'resolution', 600);
%                 print(gcf, [userpath, '/Output/Figures/', figName, '.pdf'], '-image', '-dpdf', '-r600')
                savefig(gcf, [userpath, '/Output/Figures/', figName, '.fig']);
            else
                planeNo = planeNo + 1;
            end
            
        case 'XY'
            % Initialise Figure
            if planeNo == 1
                fig = fig + 1;
%                 set(figure(fig), 'color', [1, 1, 1], 'outerPosition', [25, 25, 850, 850], 'name', figName);
                set(figure(fig), 'name', figName, 'color', [1, 1, 1], 'paperPositionMode', 'manual', 'paperUnits', 'inches', ...
                                 'paperSize', [3.45, 3.45], 'paperPosition', [0.05, 0.05, 3.35, 3.35]);
                set(gca, 'dataAspectRatio', [1, 1, 1], 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
                         'fontSize', 16, 'layer', 'top');
                lighting gouraud;
                colormap(cMap);
                hold on;
                
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
                    
                end

            end
            
            % Plot Contour
            surf(squeeze(x(:,:,2)), squeeze(y(:,:,2)), squeeze(z(:,:,2)), squeeze(scalar(:,:,2)), ...
                 'lineStyle', 'none', 'faceLighting', 'none');

            if nPlanes == 1 && ~isempty(geometry) && zLimsData > zDims(2)
                
                for i = 1:height(parts)
                    geometry.(parts{i}).boundaries.XY(:,3) = zLimsData;
                    
                    plot3(geometry.(parts{i}).boundaries.XY(:,1), ...
                          geometry.(parts{i}).boundaries.XY(:,2), ...
                          geometry.(parts{i}).boundaries.XY(:,3), ...
                          'color', 'w', 'lineStyle', '-', 'lineWidth', 1.5);
                end
                
            end
            
            % Plot Contour Lines
            if ~isempty(contourlines)
                % Convert From 'ndgrid' to 'meshgrid' Format
                x = permute(x, [2,1,3]);
                y = permute(y, [2,1,3]);
                z = permute(z, [2,1,3]);
                scalar = permute(scalar, [2,1,3]);
                
                contours = contourslice(x, y, z, scalar, [], [], zLimsData, contourlines);
                set(contours, 'edgeColor', 'w', 'lineStyle', '-', 'lineWidth', 1.5);
            end
            
            % Plot Reference Point
            if ~isempty(refPoint)
                scatter3(refPoint(1), refPoint(2), refPoint(3), 200, 'w', 'x', 'lineWidth', 2);
            end
            
            % Format Figure
            if planeNo == nPlanes

                if nPlanes == 1
                    title(figTitle, 'color', ([254, 254, 254] / 255));
                    subtitle(figSubtitle);
                    lightangle(0, 45);
                    axis on;
                    box on;
                    grid off;
                    caxis(cLims);
                    view([0, 90]);
                    xlim([xLimsPlot(1), xLimsPlot(2)]);
                    ylim([yLimsPlot(1), yLimsPlot(2)]);
                    zlim([zLimsPlot(1), zLimsPlot(2)]);
                    tickData = round((xLimsPlot(1):((xLimsPlot(2) - xLimsPlot(1)) / 5):xLimsPlot(2)), 2);
                    xticks(tickData(2:(end - 1)));
                    tickData = round((yLimsPlot(1):((yLimsPlot(2) - yLimsPlot(1)) / 5):yLimsPlot(2)), 2);
                    yticks(tickData(2:(end - 1)));
                    tickData = [];
                    zticks(tickData);
                    xtickformat('%+.2f');
                    ytickformat('%+.2f');
                    ztickformat('%+.2f');
                    
                    if normalise
                        xlabel('{\bf{x}}_{{\it{l}}}', 'fontName', 'LM Roman 12');
                        ylabel('{\bf{y}}_{{\it{l}}}', 'fontName', 'LM Roman 12');
                    else
                        xlabel('{\bf{x}}_{{\it{m}}}', 'fontName', 'LM Roman 12');
                        ylabel('{\bf{y}}_{{\it{m}}}', 'fontName', 'LM Roman 12');
                    end

                else
                    title(figTitle, 'color', ([254, 254, 254] / 255));
                    subtitle(figSubtitle);
                    lightangle(0, 45);
                    axis on;
                    box on;
                    grid off;
                    view([30, 30]);
                    xlim([xLimsPlot(1), xLimsPlot(2)]);
                    ylim([yLimsPlot(1), yLimsPlot(2)]);
                    zlim([zLimsPlot(1), zLimsPlot(2)]);
                    tickData = [];
                    xticks(tickData);
                    tickData = [];
                    yticks(tickData);
                    tickData = [];
                    zticks(tickData);
                end

%                 set(gca, 'outerPosition', [0.05, 0.05, 0.9, 0.9]);
                hold off;
                
                % Save Figure
                pause(2);
                
                exportgraphics(gcf, [userpath, '/Output/Figures/', figName, '.png'], 'resolution', 600);
%                 print(gcf, [userpath, '/Output/Figures/', figName, '.pdf'], '-image', '-dpdf', '-r600')
                savefig(gcf, [userpath, '/Output/Figures/', figName, '.fig']);
            else
                planeNo = planeNo + 1;
            end
    
    end
    
end