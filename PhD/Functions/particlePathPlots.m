%% Particle Path Plotter v2.0
% ----
% Plots Previously Processed Particle Paths
% ----
% Usage: fig = particlePathPlots();


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function fig = particlePathPlots(trackingData, fig, figName, geometry, cMap, colourVar, ...
                                 minColourVar, maxColourVar, figTitle, figSubtitle, ...
                                 xLimsPlot, yLimsPlot, zLimsPlot)
    
    % Remove Unnecessary Time Instances 
    for i = 1:height(trackingData.ID)
        index = find(isnan(trackingData.age{i}) == false);
        
        trackingData.path{i} = trackingData.path{i}(index,:);
        trackingData.age{i} = trackingData.age{i}(index);
    end
    
    % Figure Setup
    fig = fig + 1;
    set(figure(fig), 'color', [1, 1, 1], 'outerPosition', [25, 25, 850, 850], 'name', figName);
    set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
             'fontSize', 20, 'layer', 'top');
    lighting gouraud;
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
    
    for i = 1:height(trackingData.ID)
        
        if strcmp(colourVar, 'Diameter')
            pColourVar = trackingData.d(i) * 1e6;
        elseif strcmp(colourVar, 'Age')
            pColourVar = trackingData.age{i}(end);
        end
        
        cIndex = ceil(1 + (height(cMap) - 1) * ((pColourVar - minColourVar) / (maxColourVar - minColourVar)));
        
        plot3(trackingData.path{i}(:,1), trackingData.path{i}(:,2), trackingData.path{i}(:,3), ...
              'color', cMap(cIndex,:), 'lineWidth', 1.5);
        scatter3(trackingData.path{i}(end,1), trackingData.path{i}(end,2), trackingData.path{i}(end,3), ...
                 30, cMap(cIndex,:), 'filled');
    end
    
    % Figure Formatting
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
    set(gca, 'outerPosition', [0.05, 0.05, 0.9, 0.9]);
    hold off;
    
    pause(2);
    exportgraphics(gca, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);
    
end

