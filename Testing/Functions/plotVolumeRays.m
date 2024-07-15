%% Volume Ray Plotter v1.0
% ----
% Moo
% ----
% Usage: fig = plotVolumeRays();


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function fig = plotVolumeRays(positionData, fig, figName, geometry, originPoint, rayID, ...
                              figTitle, figSubtitle, xLimsPlot, yLimsPlot, zLimsPlot, figSave)
    
    % Define Target Plane Boundaries
    plane = [
             positionData(1,1), min(positionData(:,2)), min(positionData(:,3));
             positionData(1,1), min(positionData(:,2)), max(positionData(:,3));
             positionData(1,1), max(positionData(:,2)), max(positionData(:,3));
             positionData(1,1), max(positionData(:,2)), min(positionData(:,3));
             positionData(1,1), min(positionData(:,2)), min(positionData(:,3));
            ];
    
    % Initialise Figure
    fig = fig + 1;
    set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                     'outerPosition', [25, 25, 650, 650], 'units', 'pixels')
    set(gca, 'positionConstraint', 'outerPosition', 'dataAspectRatio', [1, 1, 1], ...
             'lineWidth', 2, 'fontName', 'LM Mono 12', 'fontSize', 16, 'layer', 'top');
    lighting gouraud;
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
    
    % Plot Plane Boundaries
    plot3(plane(:,1), plane(:,2), plane(:,3), 'color', ([230, 0, 126] / 255), 'lineWidth', 2);
    
    % Plot Rays
    scatter3(originPoint(1), originPoint(2), originPoint(3), 75, ([230, 0, 126] / 255), 'filled');
    
    for i = 1:height(rayID)
        plot3([originPoint(1); positionData(rayID(i),1)], ...
              [originPoint(2); positionData(rayID(i),2)], ...
              [originPoint(3); positionData(rayID(i),3)], ...
              'color', ([74, 24, 99, 127.5] / 255), 'lineWidth', 2);
    end
    
    % Format Figure
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
    hold off;
    
    % Save Figure
    pause(2);
    exportgraphics(gcf, [userpath, '/Output/Figures/', figName, '.png'], 'resolution', 600);

    if figSave
        savefig(gcf, [userpath, '/Output/Figures/', figName, '.fig']);
    end

end