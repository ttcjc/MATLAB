%% Velocity Field Plotter v1.0
% ---
% Plots Previously Processed Velocity Fields
% Usage: fig = velocityPlots(planeOrientation, planePosition, ...
%                            xLimsData, yLimsData, zLimsData, ...
%                            positionData, velocityData, ...
%                            nComponents, component, ...
%                            fig, figName, cMap, geometry, ...
%                            streamlines, modelOutline, figTime, cLims, ...
%                            xLimsPlot, yLimsPlot, zLimsPlot);
%        'planeOrientation' -> 'X', 'Y', 'Z'
%        'planePosition' -> Position in 'planeOrientation' Direction
%        '*LimsData' -> Contour Plot Limits
%        'positionData' -> Cartesian Positions
%        'velocityData' -> Three-Components of Velocity @ 'positionData'
%        'nComponents' -> Number of Components Included in Calculation
%        'components' -> Which to Include, 'u', 'v', 'w'
%        'fig' -> Figure Number
%        'figName' -> Figure Name
%        'cMap' -> Colour Map
%        'geometry' -> STL to Include in Plot
%        'streamlines' -> Include Streamlines? true/false
%        'modelOutline' -> Cartesian Positions of Model Vertices
%        'figTitle' -> Figure Title
%        'cLims' -> Colour Bounds
%        '*LimsPlot' -> Figure Limits
% ---


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function fig = velocityPlots(planeOrientation, planePosition, ...
                             xLimsData, yLimsData, zLimsData, ...
                             positionData, velocityData, ...
                             nComponents, component, ...
                             fig, figName, cMap, geometry, ...
                             streamlines, modelOutline, figTitle, cLims, ...
                             xLimsPlot, yLimsPlot, zLimsPlot)
                         
    cellSize = 0.001; % l
                         
    switch planeOrientation
        
        case 'X'
            % Generate Gridded Data            
            cellSizeX = cellSize;
            cellSizeY = (yLimsData(2) - yLimsData(1)) / round((yLimsData(2) - yLimsData(1)) / cellSize);
            cellSizeZ = (zLimsData(2) - zLimsData(1)) / round((zLimsData(2) - zLimsData(1)) / cellSize);
            
            [x, y, z] = meshgrid((planePosition - cellSizeX):cellSizeX:(planePosition + cellSizeX), ...
                                 yLimsData(1):cellSizeY:yLimsData(2), ...
                                 zLimsData(1):cellSizeZ:zLimsData(2));

            interp = scatteredInterpolant(positionData(:,2), ...
                                          positionData(:,3), ...
                                          velocityData(:,1), ...
                                          'linear', 'linear');

            u = zeros(size(x));
            u(:,2,:) = interp(y(:,2,:), z(:,2,:));

            interp = scatteredInterpolant(positionData(:,2), ...
                                          positionData(:,3), ...
                                          velocityData(:,2), ...
                                          'linear', 'linear');

            v = zeros(size(x));
            v(:,2,:) = interp(y(:,2,:), z(:,2,:));

            interp = scatteredInterpolant(positionData(:,2), ...
                                          positionData(:,3), ...
                                          velocityData(:,3), ...
                                          'linear', 'linear');

            w = zeros(size(x));
            w(:,2,:) = interp(y(:,2,:), z(:,2,:));
            
            if nComponents == 1
                velocity = eval(component);
                velocity = velocity(:,2,:);
            elseif nComponents == 2
                velocity = sqrt(v(:,2,:).^2 + w(:,2,:).^2);
            elseif nComponents == 3
                velocity = sqrt(u(:,2,:).^2 + v(:,2,:).^2 + w(:,2,:).^2);
            else
                error('Invalid Number of Velocity Components');
            end
            
            % Figure Setup
            fig = fig + 1;
            set(figure(fig), 'outerPosition', [25, 25, 850, 850], 'name', figName);
            set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
                     'fontSize', 18, 'layer', 'top');
            colormap(cMap);
            hold on;

            % Plot
            part = fieldnames(geometry);
            for j = 1:height(part)
                patch(geometry.(part{j,1}), 'faceColor', [0.5, 0.5, 0.5], ...
                                            'edgeColor', [0.5, 0.5, 0.5], ...
                                            'lineStyle', 'none');
            end

            surf(squeeze(x(:,2,:)), squeeze(y(:,2,:)), squeeze(z(:,2,:)), squeeze(velocity), ...
                 'lineStyle', 'none', 'faceLighting', 'none');
             
            if streamlines
                streams = streamslice(x, y, z, zeros(size(x)), v, w, planePosition, [], [], 'noArrows');
                set(streams, 'color', 'k');
            end
            
            if ~isempty(modelOutline)
                plot3(modelOutline(:,1), modelOutline(:,2), modelOutline(:,3), ...
                      'color', 'w', 'lineStyle', ':', 'lineWidth', 1.25);
            end
            
            % Figure Formatting
            title(figTitle);
            lightangle(90, 45);
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
            xtickformat('%+.3f');
            ytickformat('%+.3f');
            ztickformat('%+.3f');
            xT = xlabel([]);
            yT = ylabel('y_{\it{l}}');
            zT = zlabel('z_{\it{l}}');
            xT.FontName = 'LM Roman 12';
            yT.FontName = 'LM Roman 12';
            zT.FontName = 'LM Roman 12';
            hold off;

            pause(2);
            exportgraphics(gca, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);
%             savefig(fig, ['~/MATLAB/Output/Figures/', figName]);
            
        case 'Y'
            % Generate Gridded Data
            cellSizeX = (xLimsData(2) - xLimsData(1)) / round((xLimsData(2) - xLimsData(1)) / cellSize);
            cellSizeY = cellSize;
            cellSizeZ = (zLimsData(2) - zLimsData(1)) / round((zLimsData(2) - zLimsData(1)) / cellSize);
            
            [x, y, z] = meshgrid(xLimsData(1):cellSizeX:xLimsData(2), ...
                                 (planePosition - cellSizeY):cellSizeY:(planePosition + cellSizeY), ...
                                 zLimsData(1):cellSizeZ:zLimsData(2));

            interp = scatteredInterpolant(positionData(:,1), ...
                                          positionData(:,3), ...
                                          velocityData(:,1), ...
                                          'linear', 'linear');

            u = zeros(size(x));
            u(2,:,:) = interp(x(2,:,:), z(2,:,:));

            interp = scatteredInterpolant(positionData(:,1), ...
                                          positionData(:,3), ...
                                          velocityData(:,2), ...
                                          'linear', 'linear');

            v = zeros(size(x));
            v(2,:,:) = interp(x(2,:,:), z(2,:,:));

            interp = scatteredInterpolant(positionData(:,1), ...
                                          positionData(:,3), ...
                                          velocityData(:,3), ...
                                          'linear', 'linear');

            w = zeros(size(x));
            w(2,:,:) = interp(x(2,:,:), z(2,:,:));
            
            if nComponents == 1
                velocity = eval(component);
                velocity = velocity(2,:,:);
            elseif nComponents == 2
                velocity = sqrt(u(2,:,:).^2 + w(2,:,:).^2);
            elseif nComponents == 3
                velocity = sqrt(u(2,:,:).^2 + v(2,:,:).^2 + w(2,:,:).^2);
            else
                error('Invalid Number of Velocity Components');
            end
            
            % Figure Setup
            fig = fig + 1;
            set(figure(fig), 'outerPosition', [25, 25, 850, 850], 'name', figName);
            set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
                     'fontSize', 18, 'layer', 'top');
            colormap(cMap);
            hold on;

            % Plot
            part = fieldnames(geometry);
            for j = 1:height(part)
                patch(geometry.(part{j,1}), 'faceColor', [0.5, 0.5, 0.5], ...
                                            'edgeColor', [0.5, 0.5, 0.5], ...
                                            'lineStyle', 'none');
            end

            surf(squeeze(x(2,:,:)), squeeze(y(2,:,:)), squeeze(z(2,:,:)), squeeze(velocity), ...
                 'lineStyle', 'none', 'faceLighting', 'none');
             
            if streamlines
                streams = streamslice(x, y, z, u, v, w, [], positionData(1,2), [], ...
                                      1.5, 'noArrows', 'linear');
                set(streams, 'color', 'k');
            end
            
            % Figure Formatting
            title(figTitle);
            lightangle(0, 45);
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
            xtickformat('%+.3f');
            ytickformat('%+.3f');
            ztickformat('%+.3f');
            xT = xlabel('x_{\it{l}}');
            yT = ylabel([]);
            zT = zlabel('z_{\it{l}}');
            xT.FontName = 'LM Roman 12';
            yT.FontName = 'LM Roman 12';
            zT.FontName = 'LM Roman 12';
            hold off;

            pause(2);
            exportgraphics(gca, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);
%             savefig(fig, ['~/MATLAB/Output/Figures/', figName]);
            
        case 'Z'
            % Generate Gridded Data
            cellSizeX = (xLimsData(2) - xLimsData(1)) / round((xLimsData(2) - xLimsData(1)) / cellSize);
            cellSizeY = (yLimsData(2) - yLimsData(1)) / round((yLimsData(2) - yLimsData(1)) / cellSize);
            cellSizeZ = cellSize;
            
            [x, y, z] = meshgrid(xLimsData(1):cellSizeX:xLimsData(2), ...
                                 yLimsData(1):cellSizeY:yLimsData(2), ...
                                 (planePosition - cellSizeZ):cellSizeZ:(planePosition + cellSizeZ));

            interp = scatteredInterpolant(positionData(:,1), ...
                                          positionData(:,2), ...
                                          velocityData(:,1), ...
                                          'linear', 'linear');

            u = zeros(size(x));
            u(:,:,2) = interp(x(:,:,2), y(:,:,2));

            interp = scatteredInterpolant(positionData(:,1), ...
                                          positionData(:,2), ...
                                          velocityData(:,2), ...
                                          'linear', 'linear');

            v = zeros(size(x));
            v(:,:,2) = interp(x(:,:,2), y(:,:,2));

            interp = scatteredInterpolant(positionData(:,1), ...
                                          positionData(:,2), ...
                                          velocityData(:,3), ...
                                          'linear', 'linear');

            w = zeros(size(x));
            w(:,:,2) = interp(x(:,:,2), y(:,:,2));
            
            if nComponents == 1
                velocity = eval(component);
                velocity = velocity(:,:,2);
            elseif nComponents == 2
                velocity = sqrt(u(:,:,2).^2 + v(:,:,2).^2);
            elseif nComponents == 3
                velocity = sqrt(u(:,:,2).^2 + v(:,:,2).^2 + w(:,:,2).^2);
            else
                error('Invalid Number of Velocity Components');
            end
            
            % Figure Setup
            fig = fig + 1;
            set(figure(fig), 'outerPosition', [25, 25, 850, 850], 'name', figName);
            set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
                     'fontSize', 18, 'layer', 'top');
            colormap(cMap);
            hold on;

            % Plot
            part = fieldnames(geometry);
            for j = 1:height(part)
                patch(geometry.(part{j,1}), 'faceColor', [0.5, 0.5, 0.5], ...
                                            'edgeColor', [0.5, 0.5, 0.5], ...
                                            'lineStyle', 'none');
            end

            surf(squeeze(x(:,:,2)), squeeze(y(:,:,2)), squeeze(z(:,:,2)), squeeze(velocity), ...
                 'lineStyle', 'none', 'faceLighting', 'none');
             
            if streamlines
                streams = streamslice(x, y, z, u, v, zeros(size(x)), [], [], planePosition, 'noArrows');
                set(streams, 'color', 'k');
            end
            
            % Figure Formatting
            title(figTitle);
            lightangle(0, 45);
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
            xtickformat('%+.3f');
            ytickformat('%+.3f');
            ztickformat('%+.3f');
            xT = xlabel('x_{\it{l}}');
            yT = ylabel('y_{\it{l}}');
            zT = zlabel([]);
            xT.FontName = 'LM Roman 12';
            yT.FontName = 'LM Roman 12';
            zT.FontName = 'LM Roman 12';
            hold off;

            pause(2);
            exportgraphics(gca, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);
%             savefig(fig, ['~/MATLAB/Output/Figures/', figName]);
            
    end
    
end