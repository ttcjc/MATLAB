%% Particle Path Plotter v2.2
% ----
% Plots Previously Processed Particle Paths
% ----
% Usage: fig = plotParticlePaths(trackingData, fig, figName, geometry, cMap, colourVar, ...
%                                minColourVar, maxColourVar, figTitle, viewAngle, ...
%                                xLimsPlot, yLimsPlot, zLimsPlot, figSave);
%
%        'trackingData' -> Discrete Particle Paths in Cartesian Form
%        'fig'          -> Figure Number
%        'figName'      -> Figure Name
%        'geometry'     -> STL(s) to Include in Plot
%        'cMap'         -> Colour Map
%        'colourVar'    -> Particle Property Used to Colour Tracks
%        'minColourVar' -> Minimum Possible Value of 'colourVar'
%        'maxColourVar' -> Maximum Possible Value of 'colourVar'
%        'figTitle'     -> Figure Title
%        'viewAngle'    -> Default Viewing Angle
%        '*LimsPlot'    -> 3D Axes Limits
%        'figSave'      -> Save .fig File [True/False]


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats
% v2.1 - Rename and Minor Formatting Updates
% v2.2 - Update To Ensure Consistent Figure Sizes


%% Main Function

function fig = plotParticlePaths(trackingData, fig, figName, geometry, cMap, colourVar, ...
                                 minColourVar, maxColourVar, figTitle, viewAngle, ...
                                 xLimsPlot, yLimsPlot, zLimsPlot, figSave)
    
    % Initialise Figure
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
    
    % Plot Particle Tracks
    for i = 1:height(trackingData.ID)
        
        if strcmp(colourVar, 'Diameter')
            pColourVar = trackingData.d(i);
        elseif strcmp(colourVar, 'Age')
            pColourVar = trackingData.age{i}(end);
        end
        
        cIndex = ceil(1 + (height(cMap) - 1) * ((pColourVar - minColourVar) / (maxColourVar - minColourVar)));
        
        plot3(trackingData.path{i}(:,1), trackingData.path{i}(:,2), trackingData.path{i}(:,3), ...
              'color', cMap(cIndex,:), 'lineWidth', 2);
    end
    clear i;
    
    % Format Figure
    title('{-----}', 'interpreter', 'latex');
    subtitle(figTitle);
    lightangle(0, 45);
    axis on;
    box on;
    grid off
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
    tightInset = get(gca, 'TightInset');
    set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                               (tightInset(2) + 0.00625), ...
                               (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                               (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
    pause(0.5);
    hold off;
    
    % Save Figure
    print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');
    
    if figSave
        savefig(gcf, [userpath, '/Output/Figures/', figName, '.fig']);
    end
    
end