%% Vector Field Plotter v1.0
% ----
% Plots Previously Processed Vector Fields
% ----
% Usage: fig = vectorPlots(xLimsPlot, yLimsPlot, zLimsPlot, xLimsData, yLimsData, zLimsData, ...
%                          planeOrientation, planePosition, positionData, vectorData, ...
%                          nComponents, component, fig, figName, cMap, geometry, streamlines, ...
%                          xDims, yDims, zDims, figTitle, figSubtitle, cLims, normalise);
%        'planeOrientation' -> ['X', 'Y', 'Z']
%        'caseType'         -> Case Format for Setting Figure Limits
%        'normalise'        -> Normalise Dimensions [True/False]
%        'precision'        -> Rounding Precision
%        '*Lims'            -> Contour Plot Limits
%        'planePosition'    -> Cartesian Position in 'planeOrientation' Direction
%        'positionData'     -> Cartesian Positions of Data Points
%        'vectorData'       -> Three-Components of Vector @ 'positionData' Points
%        'nComponents'      -> Number of Components Included in Contour
%        'component'        -> Which Component(s) to Include, ['u', 'v', 'w']
%        'geometry'         -> STL to Include in Plot
%        'fig'              -> Figure Number
%        'figName'          -> Figure Name
%        'cMap'             -> Colour Map
%        'streamlines'      -> Include Streamlines [True/False]
%        '*Dims'            -> Simple Bounding Box of Geometry
%        'figTitle'         -> Leave Blank ('-') for Formatting Purposes
%        'figSubtitle'      -> Figure Title
%        'cLims'            -> Colour Map Limits


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function fig = vectorPlots(xLimsPlot, yLimsPlot, zLimsPlot, xLimsData, yLimsData, zLimsData, ...
                           planeOrientation, planePosition, positionData, vectorData, ...
                           nComponents, component, fig, figName, cMap, geometry, streamlines, ...
                           xDims, yDims, zDims, figTitle, figSubtitle, cLims, normalise)
    
    cellSize = 1e-3; % [m or l]
        
    % Set Plot Limits
    xLimsPlot(1) = floor(xLimsPlot(1) / 0.1) * 0.1;
    xLimsPlot(2) = floor(xLimsPlot(2) / 0.1) * 0.1;

    yLimsPlot(1) = ceil(yLimsPlot(1) / 0.1) * 0.1;
    yLimsPlot(2) = floor(yLimsPlot(2) / 0.1) * 0.1;

    zLimsPlot(1) = ceil(zLimsPlot(1) / 0.1) * 0.1;
    zLimsPlot(2) = floor(zLimsPlot(2) / 0.1) * 0.1;
    
    % Format Data
    switch planeOrientation
        
        case 'X'
            % Generate Gridded Data
            cellSizeX = cellSize;
            cellSizeY = (yLimsData(2) - yLimsData(1)) / round((yLimsData(2) - yLimsData(1)) / cellSize);
            cellSizeZ = (zLimsData(2) - zLimsData(1)) / round((zLimsData(2) - zLimsData(1)) / cellSize);
            
            [x, y, z] = meshgrid((planePosition - cellSizeX):cellSizeX:(planePosition + cellSizeX), ...
                                 yLimsData(1):cellSizeY:yLimsData(2), ...
                                 zLimsData(1):cellSizeZ:zLimsData(2));
            
            interp = scatteredInterpolant(positionData(:,2), positionData(:,3), vectorData(:,1), ...
                                          'linear', 'linear');
            
            u = zeros(size(x));
            u(:,2,:) = interp(y(:,2,:), z(:,2,:));
            
            interp = scatteredInterpolant(positionData(:,2), positionData(:,3), vectorData(:,2), ...
                                          'linear', 'linear');
            
            v = zeros(size(x));
            v(:,2,:) = interp(y(:,2,:), z(:,2,:));
            
            interp = scatteredInterpolant(positionData(:,2), positionData(:,3), vectorData(:,3), ...
                                          'linear', 'linear');
            
            w = zeros(size(x));
            w(:,2,:) = interp(y(:,2,:), z(:,2,:));
            
            % Evaluate Vector Components for Contour
            if nComponents == 1
                vector = eval(component);
                vector = vector(:,2,:);
            elseif nComponents == 2
                vector = sqrt(v(:,2,:).^2 + w(:,2,:).^2);
            elseif nComponents == 3
                vector = sqrt(u(:,2,:).^2 + v(:,2,:).^2 + w(:,2,:).^2);
            else
                error('Invalid Number of Vector Components');
            end
            
        case 'Y'
            % Generate Gridded Data
            cellSizeX = (xLimsData(2) - xLimsData(1)) / round((xLimsData(2) - xLimsData(1)) / cellSize);
            cellSizeY = cellSize;
            cellSizeZ = (zLimsData(2) - zLimsData(1)) / round((zLimsData(2) - zLimsData(1)) / cellSize);
            
            [x, y, z] = meshgrid(xLimsData(1):cellSizeX:xLimsData(2), ...
                                 (planePosition - cellSizeY):cellSizeY:(planePosition + cellSizeY), ...
                                 zLimsData(1):cellSizeZ:zLimsData(2));
            
            interp = scatteredInterpolant(positionData(:,1), positionData(:,3), vectorData(:,1), ...
                                          'linear', 'linear');
            
            u = zeros(size(x));
            u(2,:,:) = interp(x(2,:,:), z(2,:,:));
            
            interp = scatteredInterpolant(positionData(:,1), positionData(:,3), vectorData(:,2), ...
                                          'linear', 'linear');
            
            v = zeros(size(x));
            v(2,:,:) = interp(x(2,:,:), z(2,:,:));
            
            interp = scatteredInterpolant(positionData(:,1), positionData(:,3), vectorData(:,3), ...
                                          'linear', 'linear');
            
            w = zeros(size(x));
            w(2,:,:) = interp(x(2,:,:), z(2,:,:));
            
            % Evaluate Vector Components for Contour
            if nComponents == 1
                vector = eval(component);
                vector = vector(2,:,:);
            elseif nComponents == 2
                vector = sqrt(u(2,:,:).^2 + w(2,:,:).^2);
            elseif nComponents == 3
                vector = sqrt(u(2,:,:).^2 + v(2,:,:).^2 + w(2,:,:).^2);
            else
                error('Invalid Number of Vector Components');
            end
            
        case 'Z'
            % Generate Gridded Data
            cellSizeX = (xLimsData(2) - xLimsData(1)) / round((xLimsData(2) - xLimsData(1)) / cellSize);
            cellSizeY = (yLimsData(2) - yLimsData(1)) / round((yLimsData(2) - yLimsData(1)) / cellSize);
            cellSizeZ = cellSize;
            
            [x, y, z] = meshgrid(xLimsData(1):cellSizeX:xLimsData(2), ...
                                 yLimsData(1):cellSizeY:yLimsData(2), ...
                                 (planePosition - cellSizeZ):cellSizeZ:(planePosition + cellSizeZ));
            
            interp = scatteredInterpolant(positionData(:,1), positionData(:,2), vectorData(:,1), ...
                                          'linear', 'linear');
            
            u = zeros(size(x));
            u(:,:,2) = interp(x(:,:,2), y(:,:,2));
            
            interp = scatteredInterpolant(positionData(:,1), positionData(:,2), vectorData(:,2), ...
                                          'linear', 'linear');
            
            v = zeros(size(x));
            v(:,:,2) = interp(x(:,:,2), y(:,:,2));
            
            interp = scatteredInterpolant(positionData(:,1), positionData(:,2), vectorData(:,3), ...
                                          'linear', 'linear');
            
            w = zeros(size(x));
            w(:,:,2) = interp(x(:,:,2), y(:,:,2));
            
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
    switch planeOrientation
        
        case 'X'
            % Figure Setup
            fig = fig + 1;
            set(figure(fig), 'color', [1, 1, 1], 'outerPosition', [25, 25, 850, 850], 'name', figName);
            set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
                     'fontSize', 20, 'layer', 'top');
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
            
            surf(squeeze(x(:,2,:)), squeeze(y(:,2,:)), squeeze(z(:,2,:)), squeeze(vector), ...
                 'lineStyle', 'none', 'faceLighting', 'none');
             
             if streamlines
                streams = streamslice(x, y, z, zeros(size(x)), v, w, planePosition, [], [], ...
                                      2, 'arrows', 'linear');
                set(streams, 'color', 'k');
             end
            
            if planePosition < xDims(1) || planePosition > xDims(2)

                for i = 1:height(parts)
                    geometry.(parts{i}).boundaries.X(:,1) = planePosition;

                    plot3(geometry.(parts{i}).boundaries.X(:,1), ...
                          geometry.(parts{i}).boundaries.X(:,2), ...
                          geometry.(parts{i}).boundaries.X(:,3), ...
                          'color', ([230, 0, 126] / 255), 'lineStyle', '-', 'lineWidth', 1.5);
                end

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
            tickData = yLimsPlot(1):((yLimsPlot(2) - yLimsPlot(1)) / 5):yLimsPlot(2);
            yticks(tickData(2:end-1));
            tickData = zLimsPlot(1):((zLimsPlot(2) - zLimsPlot(1)) / 5):zLimsPlot(2);
            zticks(tickData(2:end-1));
            xtickformat('%+.2f');
            ytickformat('%+.2f');
            ztickformat('%+.2f');
            
            if normalise
                xT = xlabel([]);
                yT = ylabel('y_{\it{l}}');
                zT = zlabel('z_{\it{l}}');
            else
                xT = xlabel([]);
                yT = ylabel('y_{\it{m}}');
                zT = zlabel('z_{\it{m}}');
            end
            
            xT.FontName = 'LM Roman 12';
            yT.FontName = 'LM Roman 12';
            zT.FontName = 'LM Roman 12';
            hold off;
            
            pause(2);
            exportgraphics(gcf, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);
            
        case 'Y'
            % Figure Setup
            fig = fig + 1;
            set(figure(fig), 'color', [1, 1, 1], 'outerPosition', [25, 25, 850, 850], 'name', figName);
            set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
                     'fontSize', 20, 'layer', 'top');
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
            
            surf(squeeze(x(2,:,:)), squeeze(y(2,:,:)), squeeze(z(2,:,:)), squeeze(vector), ...
                 'lineStyle', 'none', 'faceLighting', 'none');
             
             if streamlines
                streams = streamslice(x, y, z, u, zeros(size(x)), w, [], planePosition, [], ...
                                      2, 'arrows', 'linear');
                set(streams, 'color', 'k');
             end
            
            if planePosition < yDims(1)

                for i = 1:height(parts)
                    geometry.(parts{i}).boundaries.Y(:,2) = planePosition;

                    plot3(geometry.(parts{i}).boundaries.Y(:,1), ...
                          geometry.(parts{i}).boundaries.Y(:,2), ...
                          geometry.(parts{i}).boundaries.Y(:,3), ...
                          'color', ([230, 0, 126] / 255), 'lineStyle', '-', 'lineWidth', 1.5);
                end

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
            tickData = xLimsPlot(1):((xLimsPlot(2) - xLimsPlot(1)) / 5):xLimsPlot(2);
            xticks(tickData(2:end-1));
            tickData = [];
            yticks(tickData);
            tickData = zLimsPlot(1):((zLimsPlot(2) - zLimsPlot(1)) / 5):zLimsPlot(2);
            zticks(tickData(2:end-1));
            xtickformat('%+.2f');
            ytickformat('%+.2f');
            ztickformat('%+.2f');
            
            if normalise
                xT = xlabel('x_{\it{l}}');
                yT = ylabel([]);
                zT = zlabel('z_{\it{l}}');
            else
                xT = xlabel('x_{\it{m}}');
                yT = ylabel([]);
                zT = zlabel('z_{\it{m}}');
            end
            
            xT.FontName = 'LM Roman 12';
            yT.FontName = 'LM Roman 12';
            zT.FontName = 'LM Roman 12';
            hold off;
            
            pause(2);
            exportgraphics(gcf, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);

        case 'Z'
            % Figure Setup
            fig = fig + 1;
            set(figure(fig), 'color', [1, 1, 1], 'outerPosition', [25, 25, 850, 850], 'name', figName);
            set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
                     'fontSize', 20, 'layer', 'top');
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
            
            surf(squeeze(x(:,:,2)), squeeze(y(:,:,2)), squeeze(z(:,:,2)), squeeze(vector), ...
                 'lineStyle', 'none', 'faceLighting', 'none');
             
             if streamlines
                streams = streamslice(x, y, z, u, v, zeros(size(x)), [], [], planePosition, ...
                                      2, 'arrows', 'linear');
                set(streams, 'color', 'k');
             end
            
            if planePosition > zDims(2)

                for i = 1:height(parts)
                    geometry.(parts{i}).boundaries.Z(:,3) = planePosition;

                    plot3(geometry.(parts{i}).boundaries.Z(:,1), ...
                          geometry.(parts{i}).boundaries.Z(:,2), ...
                          geometry.(parts{i}).boundaries.Z(:,3), ...
                          'color', ([230, 0, 126] / 255), 'lineStyle', '-', 'lineWidth', 1.5);
                end

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
            tickData = xLimsPlot(1):((xLimsPlot(2) - xLimsPlot(1)) / 5):xLimsPlot(2);
            xticks(tickData(2:end-1));
            tickData = yLimsPlot(1):((yLimsPlot(2) - yLimsPlot(1)) / 5):yLimsPlot(2);
            yticks(tickData(2:end-1));
            tickData = [];
            zticks(tickData);
            xtickformat('%+.2f');
            ytickformat('%+.2f');
            ztickformat('%+.2f');
            
            if normalise
                xT = xlabel('z_{\it{l}}');
                yT = ylabel('y_{\it{l}}');
                zT = zlabel([]);
            else
                xT = xlabel('z_{\it{m}}');
                yT = ylabel('y_{\it{m}}');
                zT = zlabel([]);

            end
            xT.FontName = 'LM Roman 12';
            yT.FontName = 'LM Roman 12';
            zT.FontName = 'LM Roman 12';
            hold off;
            
            pause(2);
            exportgraphics(gcf, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);
    
    end       