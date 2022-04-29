%% Vector Field Plotter v1.0
% ----
% Plots Previously Processed Vector Fields
% ----
% Usage: fig = vectorPlots();


%% Changelog

% v1.0 - Initial Commit


%% Supported Fields


%% Main Function

function fig = vectorPlots(planeOrientation, caseType, normalise, precision, ...
                           xLimsData, yLimsData, zLimsData, planePosition, positionData, vectorData, ...
                           nComponents, component, modelOutline, geometry, fig, figName, ...
                           cMap, streamlines, figTitle, cLims)
    
    cellSize = 1e-3; % [m or l]
    
    % Format Data
    switch planeOrientation
        
        case 'X'
            % Set Plot Limits and Normalise
            if contains(caseType, ["Run_Test", "Windsor", "Varney"])
                xLimsPlot = [0.31875; 3.61525]; % [m]
                yLimsPlot = [-0.4945; 0.4945];
                zLimsPlot = [0; 0.639];
                
                if normalise
                    xLimsPlot = round(xLimsPlot / 1.044, precision);
                    xLimsPlot(1) = floor(xLimsPlot(1) / 0.05) * 0.05;
                    xLimsPlot(2) = ceil(xLimsPlot(2) / 0.05) * 0.05;
                    
                    yLimsPlot = round(yLimsPlot / 1.044, precision);
                    yLimsPlot(1) = floor(yLimsPlot(1) / 0.05) * 0.05;
                    yLimsPlot(2) = ceil(yLimsPlot(2) / 0.05) * 0.05;

                    zLimsPlot = round(zLimsPlot / 1.044, precision);
                    zLimsPlot(1) = floor(zLimsPlot(1) / 0.05) * 0.05;
                    zLimsPlot(2) = ceil(zLimsPlot(2) / 0.05) * 0.05;
                end
                
            end
            
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
            
%             if modelOutline
%                 % Identify Model Outline
%                 parts = fieldnames(geometry);
%                 
%                 boundingGeometries = 1:height(parts);
%                 
%                 i = 1;
%                 while i <= height(parts)
%                     
%                     if ~any(round(geometry.(parts{i}).vertices(:,1), precision) == xLimsPlot(1), 'all')
%                         boundingGeometries(i) = [];
%                     else
%                         i = i + 1;
%                     end
%                     
%                 end
%                 
%                 geoPoints = [];
%                 
%                 for i = boundingGeometries
%                     geoPoints = vertcat(geoPoints, geometry.(parts{i}).vertices); %#ok<AGROW>
%                 end
%                 
%                 geoPoints = geoPoints((round(geoPoints(:,1), precision) == xLimsPlot(1)),:);
%                 
%                 modelOutline = boundary(geoPoints(:,2), geoPoints(:,3), 1);
%                 modelOutline = geoPoints(modelOutline,:);
%             else
%                 modelOutline = [];
%             end
            
        case 'Y'
            
        case 'Z'
            
    end
    
    % Present Plane
    switch planeOrientation
        
        case 'X'
            % Figure Setup
            fig = fig + 1;
            set(figure(fig), 'outerPosition', [25, 25, 850, 850], 'name', figName);
            set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
                     'fontSize', 20, 'layer', 'top');
            colormap(cMap);
            hold on;
            
            % Figure Plotting
            parts = fieldnames(geometry);

            for i = 1:height(parts)
                patch(geometry.(parts{i}), 'faceColor', ([128, 128, 128] / 255), ...
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
            
            if ~isempty(modelOutline)
                plot3(modelOutline(:,1), modelOutline(:,2), modelOutline(:,3), ...
                      'color', 'w', 'lineStyle', '-', 'lineWidth', 1.25);
            end
            
            % Figure Formatting
            title(' ', 'color', ([254, 254, 254] / 255));
            subtitle(figTitle);
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
            
        case 'Y'
            
        case 'Z'
            
    end       