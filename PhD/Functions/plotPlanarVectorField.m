%% Planar Vector Field Plotter v2.1
% ----
% Plots Previously Processed Planar Vector Fields
% ----
% Usage: [fig, planeNo] = plotPlanarVectorField(orientation, xLimsData, yLimsData, zLimsData, positionData, ...
%                                               vectorData, nComponents, component, mapPerim, ...
%                                               nPlanes, planeNo, fig, figName, cMap, geometry, ...
%                                               streamlines, xDims, yDims, zDims, figTitle, figSubtitle, ...
%                                               cLims, xLimsPlot, yLimsPlot, zLimsPlot, normalise, figSave)
% 
%        'orientation'  -> Plane Orientation ['YZ', 'XZ', 'XY']
%        '*LimsData'    -> Contour Plot Limits
%        'positionData' -> Cartesian Positions of Data Points
%        'vectorData'   -> Field Values @ 'positionData' Points
%        'nComponents'  -> 
%        'component'    ->
%        'mapPerim'     -> Map Perimeter Used When Plotting a Plane of Arbitrary Shape
%        'nPlanes'      -> Number of Planes in a Multi-Plane Figure
%        'planeNo'      -> Current Plane Number
%        'fig'          -> Figure Number
%        'figName'      -> Figure Name
%        'cMap'         -> Colour Map
%        'geometry'     -> STL(s) to Include in Plot
%        'streamlines'  -> Include Streamlines [True/False]
%        '*Dims'        -> Simple Bounding Box of Geometry
%        'figTitle'     -> Leave Blank ('-') for Formatting Purposes
%        'figSubtitle'  -> Figure Title
%        'cLims'        -> Colour Map Limits
%        '*LimsPlot'    -> 3D Axes Limits
%        'normalise'    -> Normalise Dimensions [True/False]
%        'figSave'     -> Save .fig File [True/False]


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Shifted to Using 'griddedInterpolant' and Added Support for Multi-Plane Figures)
% v2.1 - Rename and Minor Formatting Updates


%% Main Function

function [fig, planeNo] = plotPlanarVectorField(orientation, xLimsData, yLimsData, zLimsData, positionData, ...
                                                vectorData, nComponents, component, mapPerim, ...
                                                nPlanes, planeNo, fig, figName, cMap, geometry, ...
                                                streamlines, xDims, yDims, zDims, figTitle, figSubtitle, ...
                                                cLims, xLimsPlot, yLimsPlot, zLimsPlot, normalise, figSave)
    
    cellSize = 1e-3; % [m or l]
    
    % Format Data
    switch orientation
        
        case 'YZ'
            % Reshape Data for Improved Interpolation Performance
            gridShape = [height(unique(positionData(:,2))), ...
                         height(unique(positionData(:,3)))];

            y = reshape(positionData(:,2), gridShape);
            z = reshape(positionData(:,3), gridShape);

            u = reshape(vectorData(:,1), gridShape);
            v = reshape(vectorData(:,2), gridShape);
            w = reshape(vectorData(:,3), gridShape);

            % Perform Interpolation
            uInterp = griddedInterpolant(y, z, u, 'linear', 'none');
            vInterp = griddedInterpolant(y, z, v, 'linear', 'none');
            wInterp = griddedInterpolant(y, z, w, 'linear', 'none');

            % Generate 3D Surface
            cellSizeX = cellSize;
            cellSizeY = (yLimsData(2) - yLimsData(1)) / round((yLimsData(2) - yLimsData(1)) / cellSize);
            cellSizeZ = (zLimsData(2) - zLimsData(1)) / round((zLimsData(2) - zLimsData(1)) / cellSize);
            
            [x, y, z] = ndgrid((xLimsData - cellSizeX):cellSizeX:(xLimsData + cellSizeX), ...
                               yLimsData(1):cellSizeY:yLimsData(2), ...
                               zLimsData(1):cellSizeZ:zLimsData(2));
            
            % Map Data on to 3D Surface
            u = zeros(size(x));
            u(2,:,:) = uInterp(y(2,:,:), z(2,:,:));

            v = zeros(size(x));
            v(2,:,:) = vInterp(y(2,:,:), z(2,:,:));

            w = zeros(size(x));
            w(2,:,:) = wInterp(y(2,:,:), z(2,:,:));

            % Evaluate Vector Components for Contour
            if nComponents == 1
                vector = eval(component);
                vector = vector(2,:,:);
            elseif nComponents == 2
                vector = sqrt(v(2,:,:).^2 + w(2,:,:).^2);
            elseif nComponents == 3
                vector = sqrt(u(2,:,:).^2 + v(2,:,:).^2 + w(2,:,:).^2);
            end

            % Remove Data Outside Desired Map Perimeter
            if ~isempty(mapPerim)
                [indexIn, indexOn] = inpolygon(y, z, mapPerim(:,2), mapPerim(:,3));
                indexMap = double(or(indexIn, indexOn));

                vector = vector .* indexMap;
            end

        case 'XZ'
            % Reshape Data for Improved Interpolation Performance
            gridShape = [height(unique(positionData(:,1))), ...
                         height(unique(positionData(:,3)))];

            x = reshape(positionData(:,1), gridShape);
            z = reshape(positionData(:,3), gridShape);

            u = reshape(vectorData(:,1), gridShape);
            v = reshape(vectorData(:,2), gridShape);
            w = reshape(vectorData(:,3), gridShape);

            % Perform Interpolation
            uInterp = griddedInterpolant(x, z, u, 'linear', 'none');
            vInterp = griddedInterpolant(x, z, v, 'linear', 'none');
            wInterp = griddedInterpolant(x, z, w, 'linear', 'none');

            % Generate 3D Surface
            cellSizeX = (xLimsData(2) - xLimsData(1)) / round(((xLimsData(2) - xLimsData(1)) / cellSize));
            cellSizeY = cellSize;
            cellSizeZ = (zLimsData(2) - zLimsData(1)) / round((zLimsData(2) - zLimsData(1)) / cellSize);
            
            [x, y, z] = ndgrid(xLimsData(1):cellSizeX:xLimsData(2), ...
                               (yLimsData - cellSizeY):cellSizeY:(yLimsData + cellSizeY), ...
                               zLimsData(1):cellSizeZ:zLimsData(2));
            
            % Map Data on to 3D Surface
            u = zeros(size(x));
            u(:,2,:) = uInterp(x(:,2,:), z(:,2,:));

            v = zeros(size(x));
            v(:,2,:) = vInterp(x(:,2,:), z(:,2,:));

            w = zeros(size(x));
            w(:,2,:) = wInterp(x(:,2,:), z(:,2,:));

            % Evaluate Vector Components for Contour
            if nComponents == 1
                vector = eval(component);
                vector = vector(:,2,:);
            elseif nComponents == 2
                vector = sqrt(u(:,2,:).^2 + w(:,2,:).^2);
            elseif nComponents == 3
                vector = sqrt(u(:,2,:).^2 + v(:,2,:).^2 + w(:,2,:).^2);
            end

            % Remove Data Outside Desired Map Perimeter
            if ~isempty(mapPerim)
                [indexIn, indexOn] = inpolygon(x, z, mapPerim(:,1), mapPerim(:,3));
                indexMap = double(or(indexIn, indexOn));

                vector = vector .* indexMap;
            end
            
        case 'XY'
            % Reshape Data for Improved Interpolation Performance
            gridShape = [height(unique(positionData(:,1))), ...
                         height(unique(positionData(:,2)))];

            x = reshape(positionData(:,1), gridShape);
            y = reshape(positionData(:,2), gridShape);

            u = reshape(vectorData(:,1), gridShape);
            v = reshape(vectorData(:,2), gridShape);
            w = reshape(vectorData(:,3), gridShape);

            % Perform Interpolation
            uInterp = griddedInterpolant(x, y, u, 'linear', 'none');
            vInterp = griddedInterpolant(x, y, v, 'linear', 'none');
            wInterp = griddedInterpolant(x, y, w, 'linear', 'none');

            % Generate 3D Surface
            cellSizeX = (xLimsData(2) - xLimsData(1)) / round(((xLimsData(2) - xLimsData(1)) / cellSize));
            cellSizeY = (yLimsData(2) - yLimsData(1)) / round(((yLimsData(2) - yLimsData(1)) / cellSize));
            cellSizeZ = cellSize;
            
            [x, y, z] = ndgrid(xLimsData(1):cellSizeX:xLimsData(2), ...
                               yLimsData(1):cellSizeY:yLimsData(2), ...
                               (zLimsData - cellSizeZ):cellSizeZ:(zLimsData + cellSizeZ));
            
            % Map Data on to 3D Surface
            u = zeros(size(x));
            u(:,:,2) = uInterp(x(:,:,2), y(:,:,2));

            v = zeros(size(x));
            v(:,:,2) = vInterp(x(:,:,2), y(:,:,2));

            w = zeros(size(x));
            w(:,:,2) = wInterp(x(:,:,2), y(:,:,2));

            % Evaluate Vector Components for Contour
            if nComponents == 1
                vector = eval(component);
                vector = vector(:,:,2);
            elseif nComponents == 2
                vector = sqrt(v(:,:,2).^2 + w(:,:,2).^2);
            elseif nComponents == 3
                vector = sqrt(u(:,:,2).^2 + v(:,:,2).^2 + w(:,:,2).^2);
            end

            % Remove Data Outside Desired Map Perimeter
            if ~isempty(mapPerim)
                [indexIn, indexOn] = inpolygon(x, y, mapPerim(:,1), mapPerim(:,2));
                indexMap = double(or(indexIn, indexOn));

                vector = vector .* indexMap;
            end
    
    end
    
    % Present Plane
    switch orientation
        
        case 'YZ'
            % Initialise Figure
            if planeNo == 1
                fig = fig + 1;
                set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                                 'outerPosition', [25, 25, 650, 650], 'units', 'pixels')
                set(gca, 'positionConstraint', 'outerPosition', 'dataAspectRatio', [1, 1, 1], ...
                         'lineWidth', 2, 'fontName', 'LM Mono 12', 'fontSize', 16, 'layer', 'top');
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
                    clear i;
                    
                end

            end

            % Plot Vector Field
            surf(squeeze(x(2,:,:)), squeeze(y(2,:,:)), squeeze(z(2,:,:)), squeeze(vector), ...
                 'lineStyle', 'none', 'faceLighting', 'none', 'faceAlpha', 0.95);

%             % Plot Geometry Outline
%             if nPlanes == 1 && ~isempty(geometry) && (xLimsData < xDims(1) || xLimsData > xDims(2))
% 
%                 for i = 1:height(parts)
%                     geometry.(parts{i}).boundaries.YZ(:,1) = xLimsData;
% 
%                     plot3(geometry.(parts{i}).boundaries.YZ(:,1), ...
%                           geometry.(parts{i}).boundaries.YZ(:,2), ...
%                           geometry.(parts{i}).boundaries.YZ(:,3), ...
%                           'color', 'w', 'lineStyle', '-', 'lineWidth', 1.5);
%                 end
%                 clear i;
% 
%             end

            % Plot Streamlines
            if streamlines
                % Convert From 'ndgrid' to 'meshgrid' Format
                x = permute(x, [2,1,3]);
                y = permute(y, [2,1,3]);
                z = permute(z, [2,1,3]);
                v = permute(v, [2,1,3]);
                w = permute(w, [2,1,3]);
                
                streams = streamslice(x, y, z, zeros(size(x)), v, w, xLimsData, [], [], ...
                                      2, 'arrows', 'linear');
                set(streams, 'color', 'k', 'lineStyle', '-', 'lineWidth', 2);
            end
            
            % Format Figure
            if planeNo == nPlanes

                if nPlanes == 1
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
                    xticks(tickData(2:5));
                    tickData = round((zLimsPlot(1):((zLimsPlot(2) - zLimsPlot(1)) / 5):zLimsPlot(2)), 2);
                    zticks(tickData(2:5));
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
                
                hold off;
                
                % Save Figure
                pause(2);
                exportgraphics(gcf, [userpath, '/Output/Figures/', figName, '.png'], 'resolution', 600);
                
                if figSave
                    savefig(gcf, [userpath, '/Output/Figures/', figName, '.fig']);
                end
                
            else
                planeNo = planeNo + 1;
            end
            
        case 'XZ'
            % Initialise Figure
            if planeNo == 1
                fig = fig + 1;
                set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                                 'outerPosition', [25, 25, 650, 650], 'units', 'pixels')
                set(gca, 'positionConstraint', 'outerPosition', 'dataAspectRatio', [1, 1, 1], ...
                         'lineWidth', 2, 'fontName', 'LM Mono 12', 'fontSize', 16, 'layer', 'top');
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
                    clear i;
                    
                end

            end

            % Plot Vector Field
            surf(squeeze(x(:,2,:)), squeeze(y(:,2,:)), squeeze(z(:,2,:)), squeeze(vector), ...
                 'lineStyle', 'none', 'faceLighting', 'none', 'faceAlpha', 0.95);

%             % Plot Geometry Outline
%             if nPlanes == 1 && ~isempty(geometry) && yLimsData < yDims(1)
% 
%                 for i = 1:height(parts)
%                     geometry.(parts{i}).boundaries.XZ(:,1) = yLimsData;
% 
%                     plot3(geometry.(parts{i}).boundaries.XZ(:,1), ...
%                           geometry.(parts{i}).boundaries.XZ(:,2), ...
%                           geometry.(parts{i}).boundaries.XZ(:,3), ...
%                           'color', 'w', 'lineStyle', '-', 'lineWidth', 1.5);
%                 end
%                 clear i;
% 
%             end

            % Plot Streamlines
            if streamlines
                % Convert From 'ndgrid' to 'meshgrid' Format
                x = permute(x, [2,1,3]);
                y = permute(y, [2,1,3]);
                z = permute(z, [2,1,3]);
                u = permute(u, [2,1,3]);
                w = permute(w, [2,1,3]);
                
                streams = streamslice(x, y, z, u, zeros(size(x)), w, [], yLimsData, [], ...
                                      2, 'arrows', 'linear');
                set(streams, 'color', 'k', 'lineStyle', '-', 'lineWidth', 2);
            end
            
            % Format Figure
            if planeNo == nPlanes

                if nPlanes == 1
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
                    x(tickData(2:5));
                    tickData = [];
                    yticks(tickData);
                    tickData = round((zLimsPlot(1):((zLimsPlot(2) - zLimsPlot(1)) / 5):zLimsPlot(2)), 2);
                    zticks(tickData(2:5));
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
                
                hold off;
                
                % Save Figure
                pause(2);
                exportgraphics(gcf, [userpath, '/Output/Figures/', figName, '.png'], 'resolution', 600);
                
                if figSave
                    savefig(gcf, [userpath, '/Output/Figures/', figName, '.fig']);
                end
                
            else
                planeNo = planeNo + 1;
            end

        case 'XY'
            % Initialise Figure
            if planeNo == 1
                fig = fig + 1;
                set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                                 'outerPosition', [25, 25, 650, 650], 'units', 'pixels')
                set(gca, 'positionConstraint', 'outerPosition', 'dataAspectRatio', [1, 1, 1], ...
                         'lineWidth', 2, 'fontName', 'LM Mono 12', 'fontSize', 16, 'layer', 'top');
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
                    clear i;
                    
                end

            end

            % Plot Vector Field
            surf(squeeze(x(:,:,2)), squeeze(y(:,:,2)), squeeze(z(:,:,2)), squeeze(vector), ...
                 'lineStyle', 'none', 'faceLighting', 'none', 'faceAlpha', 0.95);

%             % Plot Geometry Outline
%             if nPlanes == 1 && ~isempty(geometry) && zLimsData > zDims(2)
% 
%                 for i = 1:height(parts)
%                     geometry.(parts{i}).boundaries.XZ(:,1) = zLimsData;
% 
%                     plot3(geometry.(parts{i}).boundaries.XY(:,1), ...
%                           geometry.(parts{i}).boundaries.XY(:,2), ...
%                           geometry.(parts{i}).boundaries.XY(:,3), ...
%                           'color', 'w', 'lineStyle', '-', 'lineWidth', 1.5);
%                 end
%                 clear i;
% 
%             end

            % Plot Streamlines
            if streamlines
                % Convert From 'ndgrid' to 'meshgrid' Format
                x = permute(x, [2,1,3]);
                y = permute(y, [2,1,3]);
                z = permute(z, [2,1,3]);
                u = permute(u, [2,1,3]);
                v = permute(v, [2,1,3]);
                
                streams = streamslice(x, y, z, u, v, zeros(size(x)), [], [], zLimsData, ...
                                      2, 'arrows', 'linear');
                set(streams, 'color', 'k', 'lineStyle', '-', 'lineWidth', 2);
            end
            
            % Format Figure
            if planeNo == nPlanes

                if nPlanes == 1
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
                    xticks(tickData(2:5));
                    tickData = round((yLimsPlot(1):((yLimsPlot(2) - yLimsPlot(1)) / 5):yLimsPlot(2)), 2);
                    yticks(tickData(2:5));
                    tickData = [];
                    zticks(tickData);
                    xtickformat('%+.2f');
                    ytickformat('%+.2f');
                    ztickformat('%+.2f');
        
                    if normalise
                        xlabel('{\bf{x}}_{{\it{l}}}', 'fontName', 'LM Roman 12');
                        zlabel('{\bf{y}}_{{\it{l}}}', 'fontName', 'LM Roman 12');
                    else
                        xlabel('{\bf{x}}_{{\it{m}}}', 'fontName', 'LM Roman 12');
                        zlabel('{\bf{y}}_{{\it{m}}}', 'fontName', 'LM Roman 12');
                    end

                else
                    title(figTitle, 'color', ([254, 254, 254] / 255));
                    subtitle(figSubtitle);
                    lightangle(0, 45);
                    axis on;
                    box on;
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
                
                hold off;
                
                % Save Figure
                pause(2);
                exportgraphics(gcf, [userpath, '/Output/Figures/', figName, '.png'], 'resolution', 600);
                
                if figSave
                    savefig(gcf, [userpath, '/Output/Figures/', figName, '.fig']);
                end
                
            else
                planeNo = planeNo + 1;
            end
    
    end       