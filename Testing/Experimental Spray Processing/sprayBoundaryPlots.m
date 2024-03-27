run preamble;

refValue = 0.0052166;

figSave = false; % Save .fig File(s);

%%%

% Select Relevant Geometry and Define Bounding Box
[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry(geoLoc);

parts = fieldnames(geometry);
for i = 1:height(parts)
    geometry.(parts{i}).vertices = geometry.(parts{i}).vertices / normLength;
end
clear i parts;

xDims = xDims / normLength;
yDims = yDims / normLength;
zDims = zDims / normLength;


%%

caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.0L_120s_15Hz_01/T0067_T120000_F15.mat';
caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.0L_120s_15Hz_02/T0067_T120000_F15.mat';
caseC = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.0L_120s_15Hz_03/T0067_T120000_F15.mat';
caseD = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.0L_120s_15Hz_04/T0067_T120000_F15.mat';
% caseE = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.0L_600s_03Hz_01/T0333_T600000_F3.mat';    

% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_1.0L_120s_15Hz_01/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_1.0L_120s_15Hz_02/T0067_T120000_F15.mat';
% caseC = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_1.0L_120s_15Hz_03/T0067_T120000_F15.mat';
% caseD = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_1.0L_120s_15Hz_04/T0067_T120000_F15.mat';
% caseE = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_1.0L_600s_03Hz_01/T0333_T600000_F3.mat';    

% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_1.0L_120s_15Hz_01/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_1.0L_120s_15Hz_02/T0067_T120000_F15.mat';
% caseC = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_1.0L_120s_15Hz_03/T0067_T120000_F15.mat';
% caseD = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_1.0L_120s_15Hz_04/T0067_T120000_F15.mat';
% caseE = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_1.0L_600s_03Hz_01/T0333_T600000_F3.mat';    

% Load Data
% R2R.mapDataA = load(caseA, 'mapData').mapData;
R2R.mapDataB = load(caseB, 'mapData').mapData;
% R2R.mapDataC = load(caseC, 'mapData').mapData;
% R2R.mapDataD = load(caseD, 'mapData').mapData;
% R2R.mapDataE = load(caseE, 'mapData').mapData;

sprayMaps = fieldnames(R2R);


%%

xLimsData = R2R.(sprayMaps{1}).positionGrid(1,1);
yLimsData = [min(R2R.(sprayMaps{1}).positionGrid(:,2)); max(R2R.(sprayMaps{1}).positionGrid(:,2))];
zLimsData = [min(R2R.(sprayMaps{1}).positionGrid(:,3)); max(R2R.(sprayMaps{1}).positionGrid(:,3))];

xLimsPlot = [0.3; 4.6257662];
yLimsPlot = [-0.5; 0.5];
zLimsPlot = [0; 0.5];

spatialRes = 0.5e-3 / normLength;
    
positionData = R2R.(sprayMaps{1}).positionGrid / normLength;

figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes

figName = 'Run_2_Run_Variation_Boundary';

contourlines = [0.02; 0.02];


%%

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

% Plot Contour Lines
for i = 1:height(sprayMaps)
    sprayMap = sprayMaps{i};
    
    scalarData = R2R.(sprayMap).density.mean / refValue;
    
    % Reshape Data for Improved Interpolation Performance
    gridShape = [height(unique(positionData(:,2))), ...
                 height(unique(positionData(:,3)))];

    y = reshape(positionData(:,2), gridShape);
    z = reshape(positionData(:,3), gridShape);

    scalar = reshape(scalarData, gridShape);

    % Perform Interpolation
    interp = griddedInterpolant(y, z, scalar, 'linear', 'none');

    % Generate 3D Surface
    cellSizeX = spatialRes;
    cellSizeY = (yLimsData(2) - yLimsData(1)) / round(((yLimsData(2) - yLimsData(1)) / spatialRes));
    cellSizeZ = (zLimsData(2) - zLimsData(1)) / round(((zLimsData(2) - zLimsData(1)) / spatialRes));

    [x, y, z] = ndgrid((xLimsData - cellSizeX):cellSizeX:(xLimsData + cellSizeX), ...
                       yLimsData(1):cellSizeY:yLimsData(2), ...
                       zLimsData(1):cellSizeZ:zLimsData(2));

    % Map Data on to 3D Surface
    scalar = zeros(size(x));
    scalar(2,:,:) = interp(y(2,:,:), z(2,:,:));
    x = permute(x, [2,1,3]);
    y = permute(y, [2,1,3]);
    z = permute(z, [2,1,3]);
    scalar = permute(scalar, [2,1,3]);

    contours = contourslice(x, y, z, scalar, xLimsData, [], [], contourlines);
    set(contours, 'edgeColor', graphColours(i), 'lineStyle', '-', 'lineWidth', 2);
end

% % Plot Sample Points for Convergence Testing
% for i = 1:height(index)
%     plot3(positionData(index(i),1), positionData(index(i),2), positionData(index(i),3), ...
%           'lineStyle', 'none', 'marker', 'o', 'markerSize', 7.5, 'markerFaceColor', graphColours(i), ...
%           'markerEdgeColor', graphColours(i));
% end
% clear i;

% Plot Dummy Data for Legend
dummyData = zeros([height(sprayMaps),1]);
dummyData(1) = plot(NaN, NaN, 'color', graphColours(1), 'lineWidth', 2);
dummyData(2) = plot(NaN, NaN, 'color', graphColours(2), 'lineWidth', 2);
dummyData(3) = plot(NaN, NaN, 'color', graphColours(3), 'lineWidth', 2);
dummyData(4) = plot(NaN, NaN, 'color', graphColours(4), 'lineWidth', 2);

title('{-----}', 'interpreter', 'latex');
subtitle(figTitle);
lightangle(90, 45);
axis on;
box on;
grid off;
view([90, 0]);
xlim([xLimsPlot(1), xLimsPlot(2)]);
ylim([yLimsPlot(1), yLimsPlot(2)]);
zlim([zLimsPlot(1), zLimsPlot(2)]);
tickData = [];
xticks(tickData);
tickData = (yLimsPlot(1):((yLimsPlot(2) - yLimsPlot(1)) / 5):yLimsPlot(2));
yticks(tickData(2:(end-1)));
tickData = (zLimsPlot(1):((zLimsPlot(2) - zLimsPlot(1)) / 5):zLimsPlot(2));
zticks(tickData(2:(end-1)));
xtickformat('%+.2g');
ytickformat('%+.2g');
ztickformat('%+.2g');
ylabel({'{$y_{\ell}$}'; '{-----}'}, 'interpreter', 'latex');
zlabel({'{-----}'; '{$z_{\ell}$}'}, 'interpreter', 'latex');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
legProps = legend(dummyData, {'Run 1', ...
                              'Run 2', ...
                              'Run 3', ...
                              'Run 4'}, ...
                             'location', 'northEast', 'orientation', 'vertical', 'interpreter', 'latex', ...
                             'fontSize', 18, 'box', 'off');
legProps.Position(2) = legProps.Position(2) - 0.2878;
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');