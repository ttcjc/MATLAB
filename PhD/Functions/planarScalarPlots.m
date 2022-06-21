%% Planar Scalar Field Plotter v1.0
% ----
% Plots Previously Processed Planar Scalar Fields
% ----
% Usage: fig = planarScalarPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
%                                mapPerim, fig, figName, cMap, geometry, xDims, yDims, zDims, ...
%                                CoM, figTitle, figSubtitle, cLims, xLimsPlot, yLimsPlot, zLimsPlot, normalise);
%        'orientation'  -> Plane Orientation ['YZ', 'XZ', 'XY']
%        '*LimsData'    -> Contour Plot Limits
%        'positionData' -> Cartesian Positions of Data Points
%        'scalarData'   -> Field Value @ 'positionData' Points
%        'mapPerim'     -> Map Perimeter Used When Plotting a Plane of Arbitrary Shape
%        'fig'          -> Figure Number
%        'figName'      -> Figure Name
%        'cMap'         -> Colour Map
%        'geometry'     -> STL(s) to Include in Plot
%        '*Dims'        -> Simple Bounding Box of Geometry
%        'CoM'          -> Field Centre of Mass
%        'figTitle'     -> Leave Blank ('-') for Formatting Purposes
%        'figSubtitle'  -> Figure Title
%        'cLims'        -> Colour Map Limits
%        '*LimsPlot'    -> 3D Axes Limits
%        'normalise'    -> Normalise Dimensions [True/False]


%% Changelog

% v1.0 - Initial Commit (Expanded Functionality of 'contaminantPlots.m')


%% Main Function

function fig = planarScalarPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                 mapPerim, fig, figName, cMap, geometry, xDims, yDims, zDims, ...
                                 CoM, figTitle, figSubtitle, cLims, xLimsPlot, yLimsPlot, zLimsPlot, normalise)
    
    cellSize = 0.5e-3; % [m or l]
    
    % Format Data
    switch orientation
        
        case 'YZ'
            % Generate Gridded Data
            cellSizeX = cellSize;
            cellSizeY = (yLimsData(2) - yLimsData(1)) / round(((yLimsData(2) - yLimsData(1)) / cellSize));
            cellSizeZ = (zLimsData(2) - zLimsData(1)) / round(((zLimsData(2) - zLimsData(1)) / cellSize));
            
            [x, y, z] = ndgrid((xLimsData - cellSizeX):cellSizeX:(xLimsData + cellSizeX), ...
                               yLimsData(1):cellSizeY:yLimsData(2), ...
                               zLimsData(1):cellSizeZ:zLimsData(2));
                             
            interp = scatteredInterpolant(positionData(:,2), positionData(:,3), scalarData, ...
                                          'linear', 'none');
    
            scalar = zeros(size(x));
            scalar(2,:,:) = interp(y(2,:,:), z(2,:,:));
            
            if ~isempty(mapPerim)
                [indexIn, indexOn] = inpolygon(y, z, mapPerim(:,2), mapPerim(:,3));
                indexMap = double(or(indexIn, indexOn));
                indexMap(indexMap == 0) = nan;
                
                scalar = scalar .* indexMap;
            end
            
        case 'XZ'
            % Generate Gridded Data
            cellSizeX = (xLimsData(2) - xLimsData(1)) / round(((xLimsData(2) - xLimsData(1)) / cellSize));
            cellSizeY = cellSize;
            cellSizeZ = (zLimsData(2) - zLimsData(1)) / round(((zLimsData(2) - zLimsData(1)) / cellSize));

            [x, y, z] = ndgrid(xLimsData(1):cellSizeX:xLimsData(2), ...
                               (yLimsData - cellSizeY):cellSizeY:(yLimsData + cellSizeY), ...
                               zLimsData(1):cellSizeZ:zLimsData(2));
                             
            interp = scatteredInterpolant(positionData(:,1), positionData(:,3), scalarData, ...
                                          'linear', 'none');
    
            scalar = zeros(size(x));
            scalar(:,2,:) = interp(y(:,2,:), z(:,2,:));
            
            if ~isempty(mapPerim)
                [indexIn, indexOn] = inpolygon(x, z, mapPerim(:,1), mapPerim(:,3));
                indexMap = double(or(indexIn, indexOn));
                indexMap(indexMap == 0) = nan;
                
                scalar = scalar .* indexMap;
            end
            
        case 'XY'
            % Generate Gridded Data
            cellSizeX = (xLimsData(2) - xLimsData(1)) / round(((xLimsData(2) - xLimsData(1)) / cellSize));
            cellSizeY = (yLimsData(2) - yLimsData(1)) / round(((yLimsData(2) - yLimsData(1)) / cellSize));
            cellSizeZ = cellSize;
            
            [x, y, z] = ndgrid(xLimsData(1):cellSizeX:xLimsData(2), ...
                               yLimsData(1):cellSizeY:yLimsData(2), ...
                               (zLimsData - cellSizeZ):cellSizeZ:(zLimsData + cellSizeZ));
                             
            interp = scatteredInterpolant(positionData(:,1), positionData(:,2), scalarData, ...
                                          'linear', 'none');
    
            scalar = zeros(size(x));
            scalar(:,:,2) = interp(y(:,:,2), z(:,:,2));
            
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
            % Figure Setup
            fig = fig + 1;
            set(figure(fig), 'color', [1, 1, 1], 'outerPosition', [25, 25, 850, 850], 'name', figName);
            set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
                     'fontSize', 20, 'layer', 'top');
            lighting gouraud;
            colormap(cMap);
            hold on;
            
            % Figure Plotting
            parts = fieldnames(geometry);
            
            for i = 1:height(parts)
                patch('faces', geometry.(parts{i}).faces, ...
                      'vertices', geometry.(parts{i}).vertices, ...
                      'faceColor', ([128, 128, 128] / 255), ...
                      'edgeColor', ([128, 128, 128] / 255), ...
                      'lineStyle', 'none');
            end
            
            surf(squeeze(x(2,:,:)), squeeze(y(2,:,:)), squeeze(z(2,:,:)), squeeze(scalar(2,:,:)), ...
                 'lineStyle', 'none', 'faceLighting', 'none');
             
            if xLimsData < xDims(1) || xLimsData > xDims(2)
                
                for i = 1:height(parts)
                    geometry.(parts{i}).boundaries.YZ(:,1) = xLimsData;
                    
                    plot3(geometry.(parts{i}).boundaries.YZ(:,1), ...
                          geometry.(parts{i}).boundaries.YZ(:,2), ...
                          geometry.(parts{i}).boundaries.YZ(:,3), ...
                          'color', 'w', 'lineStyle', '-', 'lineWidth', 1.5);
                end
                
            end
            
            if ~isempty(CoM)
                scatter3(CoM(1), CoM(2), CoM(3), 125, 'w', 'filled');
            end
            
            % Figure Formatting
            title(figTitle, 'color', ([254, 254, 254] / 255));
            subtitle(figSubtitle);
            lightangle(90, 45);
            axis on;
            box on;
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
            
            set(gca, 'outerPosition', [0.05, 0.05, 0.9, 0.9]);
            hold off;
            
            pause(2);
            exportgraphics(gcf, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);
            
        case 'XZ'
            % Figure Setup
            fig = fig + 1;
            set(figure(fig), 'color', [1, 1, 1], 'outerPosition', [25, 25, 850, 850], 'name', figName);
            set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
                     'fontSize', 20, 'layer', 'top');
            lighting gouraud;
            colormap(cMap);
            hold on;
            
            % Figure Plotting
            parts = fieldnames(geometry);
            
            for i = 1:height(parts)
                patch('faces', geometry.(parts{i}).faces, ...
                      'vertices', geometry.(parts{i}).vertices, ...
                      'faceColor', ([128, 128, 128] / 255), ...
                      'edgeColor', ([128, 128, 128] / 255), ...
                      'lineStyle', 'none');
            end
            
            surf(squeeze(x(:,2,:)), squeeze(y(:,2,:)), squeeze(z(:,2,:)), squeeze(scalar(:,2,:)), ...
                 'lineStyle', 'none', 'faceLighting', 'none');
             
            if yLimsData < yDims(1)
                
                for i = 1:height(parts)
                    geometry.(parts{i}).boundaries.XZ(:,2) = yLimsData;
                    
                    plot3(geometry.(parts{i}).boundaries.XZ(:,1), ...
                          geometry.(parts{i}).boundaries.XZ(:,2), ...
                          geometry.(parts{i}).boundaries.XZ(:,3), ...
                          'color', 'w', 'lineStyle', '-', 'lineWidth', 1.5);
                end
                
            end
            
            if ~isempty(CoM)
                scatter3(CoM(1), CoM(2), CoM(3), 125, 'w', 'filled');
            end
            
            % Figure Formatting
            title(figTitle, 'color', ([254, 254, 254] / 255));
            subtitle(figSubtitle);
            lightangle(0, 45);
            axis on;
            box on;
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

            set(gca, 'outerPosition', [0.05, 0.05, 0.9, 0.9]);
            hold off;
            
            pause(2);
            exportgraphics(gcf, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);
            
        case 'XY'
            % Figure Setup
            fig = fig + 1;
            set(figure(fig), 'color', [1, 1, 1], 'outerPosition', [25, 25, 850, 850], 'name', figName);
            set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
                     'fontSize', 20, 'layer', 'top');
            lighting gouraud;
            colormap(cMap);
            hold on;
            
            % Figure Plotting
            parts = fieldnames(geometry);
            
            for i = 1:height(parts)
                patch('faces', geometry.(parts{i}).faces, ...
                      'vertices', geometry.(parts{i}).vertices, ...
                      'faceColor', ([128, 128, 128] / 255), ...
                      'edgeColor', ([128, 128, 128] / 255), ...
                      'lineStyle', 'none');
            end
            
            surf(squeeze(x(:,:,2)), squeeze(y(:,:,2)), squeeze(z(:,:,2)), squeeze(scalar(:,:,2)), ...
                 'lineStyle', 'none', 'faceLighting', 'none');
             
            if yLimsData > zDims(2)
                
                for i = 1:height(parts)
                    geometry.(parts{i}).boundaries.XY(:,3) = zLimsData;
                    
                    plot3(geometry.(parts{i}).boundaries.XY(:,1), ...
                          geometry.(parts{i}).boundaries.XY(:,2), ...
                          geometry.(parts{i}).boundaries.XY(:,3), ...
                          'color', 'w', 'lineStyle', '-', 'lineWidth', 1.5);
                end
                
            end
            
            if ~isempty(CoM)
                scatter3(CoM(1), CoM(2), CoM(3), 125, 'w', 'filled');
            end
            
            % Figure Formatting
            title(figTitle, 'color', ([254, 254, 254] / 255));
            subtitle(figSubtitle);
            lightangle(0, 45);
            axis on;
            box on;
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

            set(gca, 'outerPosition', [0.05, 0.05, 0.9, 0.9]);
            hold off;
            
            pause(2);
            exportgraphics(gcf, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);
    
    end
    
end