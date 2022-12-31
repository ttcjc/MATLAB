%% Vector Field Plotter v1.0
% ----
% Plots Previously Processed Planar Vector Fields
% ----
% Usage: fig = planarVectorPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, ...
%                                vectorData, nComponents, component, fig, figName, cMap, geometry, ...
%                                streamlines, xDims, yDims, zDims, figTitle, figSubtitle, cLims, ...
%                                xLimsPlot, yLimsPlot, zLimsPlot, normalise);
%        'orientation'  -> Plane Orientation ['YZ', 'XZ', 'XY']
%        '*LimsData'    -> Contour Plot Limits
%        'positionData' -> Cartesian Positions of Data Points
%        'vectorData'   -> Field Values @ 'positionData' Points
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


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function fig = planarVectorPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, ...
                                 vectorData, nComponents, component, fig, figName, cMap, geometry, ...
                                 streamlines, xDims, yDims, zDims, figTitle, figSubtitle, cLims, ...
                                 xLimsPlot, yLimsPlot, zLimsPlot, normalise)
    
    cellSize = 0.5e-3; % [m or l]
    
    % Format Data
    switch orientation
        
        case 'YZ'
            % Generate Gridded Data
            cellSizeX = cellSize;
            cellSizeY = (yLimsData(2) - yLimsData(1)) / round((yLimsData(2) - yLimsData(1)) / cellSize);
            cellSizeZ = (zLimsData(2) - zLimsData(1)) / round((zLimsData(2) - zLimsData(1)) / cellSize);
            
            [x, y, z] = ndgrid((xLimsData - cellSizeX):cellSizeX:(xLimsData + cellSizeX), ...
                               yLimsData(1):cellSizeY:yLimsData(2), ...
                               zLimsData(1):cellSizeZ:zLimsData(2));
            
            interp = scatteredInterpolant(positionData(:,2), positionData(:,3), vectorData(:,1), ...
                                          'linear', 'none');
            
            u = zeros(size(x));
            u(2,:,:) = interp(y(2,:,:), z(2,:,:));
            u(isnan(u)) = 0;
            
            interp = scatteredInterpolant(positionData(:,2), positionData(:,3), vectorData(:,2), ...
                                          'linear', 'none');
            
            v = zeros(size(x));
            v(2,:,:) = interp(y(2,:,:), z(2,:,:));
            v(isnan(v)) = 0;
            
            interp = scatteredInterpolant(positionData(:,2), positionData(:,3), vectorData(:,3), ...
                                          'linear', 'none');
            
            w = zeros(size(x));
            w(2,:,:) = interp(y(2,:,:), z(2,:,:));
            w(isnan(w)) = 0;
            
            % Evaluate Vector Components for Contour
            if nComponents == 1
                vector = eval(component);
                vector = vector(2,:,:);
            elseif nComponents == 2
                vector = sqrt(v(2,:,:).^2 + w(2,:,:).^2);
            elseif nComponents == 3
                vector = sqrt(u(2,:,:).^2 + v(2,:,:).^2 + w(2,:,:).^2);
            else
                error('Invalid Number of Vector Components');
            end
            
        case 'XZ'
            % Generate Gridded Data
            cellSizeX = (xLimsData(2) - xLimsData(1)) / round((xLimsData(2) - xLimsData(1)) / cellSize);
            cellSizeY = cellSize;
            cellSizeZ = (zLimsData(2) - zLimsData(1)) / round((zLimsData(2) - zLimsData(1)) / cellSize);
            
            [x, y, z] = ndgrid(xLimsData(1):cellSizeX:xLimsData(2), ...
                               (yLimsData - cellSizeY):cellSizeY:(yLimsData + cellSizeY), ...
                               zLimsData(1):cellSizeZ:zLimsData(2));
            
            interp = scatteredInterpolant(positionData(:,1), positionData(:,3), vectorData(:,1), ...
                                          'linear', 'none');
            
            u = zeros(size(x));
            u(:,2,:) = interp(x(:,2,:), z(:,2,:));
            u(isnan(u)) = 0;
            
            interp = scatteredInterpolant(positionData(:,1), positionData(:,3), vectorData(:,2), ...
                                          'linear', 'none');
            
            v = zeros(size(x));
            v(:,2,:) = interp(x(:,2,:), z(:,2,:));
            v(isnan(v)) = 0;
            
            interp = scatteredInterpolant(positionData(:,1), positionData(:,3), vectorData(:,3), ...
                                          'linear', 'none');
            
            w = zeros(size(x));
            w(:,2,:) = interp(x(:,2,:), z(:,2,:));
            w(isnan(w)) = 0;
            
            % Evaluate Vector Components for Contour
            if nComponents == 1
                vector = eval(component);
                vector = vector(:,2,:);
            elseif nComponents == 2
                vector = sqrt(u(:,2,:).^2 + w(:,2,:).^2);
            elseif nComponents == 3
                vector = sqrt(u(:,2,:).^2 + v(:,2,:).^2 + w(:,2,:).^2);
            else
                error('Invalid Number of Vector Components');
            end
            
        case 'XY'
            % Generate Gridded Data
            cellSizeX = (xLimsData(2) - xLimsData(1)) / round((xLimsData(2) - xLimsData(1)) / cellSize);
            cellSizeY = (yLimsData(2) - yLimsData(1)) / round((yLimsData(2) - yLimsData(1)) / cellSize);
            cellSizeZ = cellSize;
            
            [x, y, z] = ndgrid(xLimsData(1):cellSizeX:xLimsData(2), ...
                               yLimsData(1):cellSizeY:yLimsData(2), ...
                               (zLimsData - cellSizeZ):cellSizeZ:(zLimsData + cellSizeZ));
            
            interp = scatteredInterpolant(positionData(:,1), positionData(:,2), vectorData(:,1), ...
                                          'linear', 'none');
            
            u = zeros(size(x));
            u(:,:,2) = interp(x(:,:,2), y(:,:,2));
            u(isnan(u)) = 0;
            
            interp = scatteredInterpolant(positionData(:,1), positionData(:,2), vectorData(:,2), ...
                                          'linear', 'none');
            
            v = zeros(size(x));
            v(:,:,2) = interp(x(:,:,2), y(:,:,2));
            v(isnan(v)) = 0;
            
            interp = scatteredInterpolant(positionData(:,1), positionData(:,2), vectorData(:,3), ...
                                          'linear', 'none');
            
            w = zeros(size(x));
            w(:,:,2) = interp(x(:,:,2), y(:,:,2));
            w(isnan(w)) = 0;
            
            % Evaluate Vector Components for Contour
            if nComponents == 1
                vector = eval(component);
                vector = vector(:,:,2);
            elseif nComponents == 2
                vector = sqrt(u(:,:,2).^2 + v(:,:,2).^2);
            elseif nComponents == 3
                vector = sqrt(u(:,:,2).^2 + v(:,:,2).^2 + w(:,:,2).^2);
            else
                error('Invalid Number of Vector Components');
            end
    
    end
    
    % Present Plane
    switch orientation
        
        case 'YZ'
            % Figure Setup
            fig = fig + 1;
            set(figure(fig), 'color', [1, 1, 1], 'outerPosition', [25, 25, 850, 850], 'name', figName);
            set(gca, 'dataAspectRatio', [1, 1, 1], 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
                     'fontSize', 20, 'layer', 'top');
            lighting gouraud;
            colormap(cMap);
            hold on;
            
            % Figure Plotting
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
            surf(squeeze(x(2,:,:)), squeeze(y(2,:,:)), squeeze(z(2,:,:)), squeeze(vector), ...
                 'lineStyle', 'none', 'faceLighting', 'none');
                        
            if ~isempty(geometry) && (xLimsData < xDims(1) || xLimsData > xDims(2))

                for i = 1:height(parts)
                    geometry.(parts{i}).boundaries.YZ(:,1) = xLimsData;

                    plot3(geometry.(parts{i}).boundaries.YZ(:,1), ...
                          geometry.(parts{i}).boundaries.YZ(:,2), ...
                          geometry.(parts{i}).boundaries.YZ(:,3), ...
                          'color', [0.5, 0.5, 0.5], 'lineStyle', '-', 'lineWidth', 1.5);
                end

            end
            
            if streamlines
                % Convert From 'ndgrid' to 'meshgrid' Format
                x = permute(x, [2,1,3]);
                y = permute(y, [2,1,3]);
                z = permute(z, [2,1,3]);
                v = permute(v, [2,1,3]);
                w = permute(w, [2,1,3]);
                
                streams = streamslice(x, y, z, zeros(size(x)), v, w, xLimsData, [], [], ...
                                      2, 'arrows', 'linear');
                set(streams, 'color', 'k', 'lineStyle', '-', 'lineWidth', 1);
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
            set(gca, 'dataAspectRatio', [1, 1, 1], 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
                     'fontSize', 20, 'layer', 'top');
            lighting gouraud;
            colormap(cMap);
            hold on;
            
            % Figure Plotting
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
            
            surf(squeeze(x(:,2,:)), squeeze(y(:,2,:)), squeeze(z(:,2,:)), squeeze(vector), ...
                 'lineStyle', 'none', 'faceLighting', 'none');
            
            if ~isempty(geometry) && yLimsData < yDims(1)

                for i = 1:height(parts)
                    geometry.(parts{i}).boundaries.XZ(:,2) = yLimsData;

                    plot3(geometry.(parts{i}).boundaries.XZ(:,1), ...
                          geometry.(parts{i}).boundaries.XZ(:,2), ...
                          geometry.(parts{i}).boundaries.XZ(:,3), ...
                          'color', [0.5, 0.5, 0.5], 'lineStyle', '-', 'lineWidth', 1.5);
                end

            end
            
            if streamlines
                % Convert From 'ndgrid' to 'meshgrid' Format
                x = permute(x, [2,1,3]);
                y = permute(y, [2,1,3]);
                z = permute(z, [2,1,3]);
                u = permute(u, [2,1,3]);
                w = permute(w, [2,1,3]);

                streams = streamslice(x, y, z, u, zeros(size(x)), w, [], yLimsData, [], ...
                                      2, 'arrows', 'linear');
                set(streams, 'color', 'k', 'lineStyle', '-', 'lineWidth', 1);
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
            set(gca, 'dataAspectRatio', [1, 1, 1], 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
                     'fontSize', 20, 'layer', 'top');
            lighting gouraud;
            colormap(cMap);
            hold on;
            
            % Figure Plotting
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
            
            surf(squeeze(x(:,:,2)), squeeze(y(:,:,2)), squeeze(z(:,:,2)), squeeze(vector), ...
                 'lineStyle', 'none', 'faceLighting', 'none');

             if ~isempty(geometry) && zLimsData > zDims(2)

                for i = 1:height(parts)
                    geometry.(parts{i}).boundaries.XY(:,3) = zLimsData;

                    plot3(geometry.(parts{i}).boundaries.XY(:,1), ...
                          geometry.(parts{i}).boundaries.XY(:,2), ...
                          geometry.(parts{i}).boundaries.XY(:,3), ...
                          'color', [0.5, 0.5, 0.5], 'lineStyle', '-', 'lineWidth', 1.5);
                end

             end
            
            if streamlines
                % Convert From 'ndgrid' to 'meshgrid' Format
                x = permute(x, [2,1,3]);
                y = permute(y, [2,1,3]);
                z = permute(z, [2,1,3]);
                u = permute(u, [2,1,3]);
                v = permute(v, [2,1,3]);
                
                streams = streamslice(x, y, z, u, v, zeros(size(x)), [], [], zLimsData, ...
                                      2, 'arrows', 'linear');
                set(streams, 'color', 'k', 'lineStyle', '-', 'lineWidth', 1);
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