%% Particle Path Plotter v2.1
% ----
% Plots Previously Processed Particle Paths
% ----
% Usage: fig = plotParticlePaths(trackingData, fig, figName, geometry, cMap, colourVar, ...
%                                minColourVar, maxColourVar, figTitle, figSubtitle, ...
%                                xLimsPlot, yLimsPlot, zLimsPlot, figSave)
% 
%        'trackingData' -> Discrete Particle Paths in Cartesian Form
%        'fig'          -> Figure Number
%        'figName'      -> Figure Name
%        'geometry'     -> STL(s) to Include in Plot
%        'cMap'         -> Colour Map
%        'colourVar'    -> Particle Property Used to Colour Tracks
%        'minColourVar' -> Minimum Possible Value of 'colourVar'
%        'maxColourVar' -> Maximum Possible Value of 'colourVar'
%        'figTitle'     -> Leave Blank ('-') for Formatting Purposes
%        'figSubtitle'  -> Figure Title
%        '*LimsPlot'    -> 3D Axes Limits
%        'figSave'      -> Save .fig File [True/False]


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats
% v2.1 - Rename and Minor Formatting Updates


%% Main Function

function fig = plotParticlePaths(trackingData, fig, figName, geometry, cMap, colourVar, ...
                                 minColourVar, maxColourVar, figTitle, figSubtitle, ...
                                 xLimsPlot, yLimsPlot, zLimsPlot, figSave)
    
    % Remove Unnecessary Time Instances 
    for i = 1:height(trackingData.ID)
        index = find(isnan(trackingData.age{i}) == false);
        
        trackingData.path{i} = trackingData.path{i}(index,:);
        trackingData.age{i} = trackingData.age{i}(index);
    end
    
    % Initialise Figure
    fig = fig + 1;
    set(figure(fig), 'name', figName, 'color', [1, 1, 1], 'paperPositionMode', 'manual', 'paperUnits', 'inches', ...
                     'paperSize', [3.45, 3.45], 'paperPosition', [0.05, 0.05, 3.35, 3.35]);
    set(gca, 'dataAspectRatio', [1, 1, 1], 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
             'fontSize', 20, 'layer', 'top');
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
    
    % Plot Particle Tracks
    for i = 1:height(trackingData.ID)
        
        if strcmp(colourVar, 'Diameter')
            pColourVar = trackingData.d(i) * 1e6;
        elseif strcmp(colourVar, 'Age')
            pColourVar = trackingData.age{i}(end);
        end
        
        cIndex = ceil(1 + (height(cMap) - 1) * ((pColourVar - minColourVar) / (maxColourVar - minColourVar)));
        
        plot3(trackingData.path{i}(:,1), trackingData.path{i}(:,2), trackingData.path{i}(:,3), ...
              'color', cMap(cIndex,:), 'lineWidth', 1.5);
    end
    clear i;
    
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