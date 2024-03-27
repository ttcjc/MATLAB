run preamble;

figSave = false; % Save .fig File(s);

%%%

% Select Relevant Geometry and Define Bounding Box
[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry(geoLoc);


%%
% 
% clc;
% close all;


%%

% 1.0L Config Repeatability
% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.0L_120s_15Hz_02/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.0L_120s_15Hz_01/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.0L_120s_15Hz_03/T0067_T120000_F15.mat';

% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_1.0L_120s_15Hz_01/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_1.0L_120s_15Hz_02/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_1.0L_120s_15Hz_03/T0067_T120000_F15.mat';

% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_1.0L_120s_15Hz_02/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_1.0L_120s_15Hz_01/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_1.0L_120s_15Hz_03/T0067_T120000_F15.mat';


% Per-Plane Config Variations
% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.0L_120s_15Hz_02/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_1.0L_120s_15Hz_01/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_1.0L_120s_15Hz_02/T0067_T120000_F15.mat';

% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.5L_120s_15Hz_01/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_1.5L_120s_15Hz_02/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_1.5L_120s_15Hz_03/T0067_T120000_F15.mat';

caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_2.0L_120s_15Hz_03/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_2.0L_120s_15Hz_03/T0067_T120000_F15.mat';
caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_2.0L_120s_15Hz_02/T0067_T120000_F15.mat';


% Other Stuff
% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.0L_120s_15Hz_02/T0067_T120000_F15.mat';
% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.5L_120s_15Hz_01/T0067_T120000_F15.mat';
% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_2.0L_120s_15Hz_03/T0067_T120000_F15.mat';

% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.0L_120s_15Hz_02/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.5L_120s_15Hz_01/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_2.0L_120s_15Hz_03/T0067_T120000_F15.mat';


% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_1.0L_120s_15Hz_01/T0067_T120000_F15.mat';
% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_1.5L_120s_15Hz_02/T0067_T120000_F15.mat';
% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_2.0L_120s_15Hz_03/T0067_T120000_F15.mat';

% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_1.0L_120s_15Hz_01/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_1.5L_120s_15Hz_02/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/ST_2.0L_120s_15Hz_03/T0067_T120000_F15.mat';


% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_1.0L_120s_15Hz_02/T0067_T120000_F15.mat';
% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_1.5L_120s_15Hz_03/T0067_T120000_F15.mat';
% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_2.0L_120s_15Hz_02/T0067_T120000_F15.mat';

% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_1.0L_120s_15Hz_02/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_1.5L_120s_15Hz_03/T0067_T120000_F15.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/RSST_2.0L_120s_15Hz_02/T0067_T120000_F15.mat';


% Load Data
mapDataA = load(caseA, 'mapData').mapData;
mapDataB = load(caseB, 'mapData').mapData;

caseID_A = load(caseA, 'caseID').caseID;
caseID_B = load(caseB, 'caseID').caseID;


%%

% correlation = corr(mapDataA.density.RMS, mapDataB.density.RMS)

% valsA = [min(mapDataA.density.mean), max(mapDataA.density.mean), mean(mapDataA.density.mean(mapDataA.density.mean > 1e-6)), mean(mapDataA.density.RMS(mapDataA.density.mean > 1e-6))]
% valsB = [min(mapDataB.density.mean), max(mapDataB.density.mean), mean(mapDataB.density.mean(mapDataB.density.mean > 1e-6)), mean(mapDataB.density.RMS(mapDataB.density.mean > 1e-6))]

% RMSerror = (sqrt((1 / height(mapDataA.time)) * sum((mapDataA.density.mean - mapDataB.density.mean).^2))) / mean(mapDataA.density.mean)


%%

mapDataA.density.mean = mapDataA.density.mean / 0.0052166;
mapDataB.density.mean = mapDataB.density.mean / 0.0052166;

meanDelta = mapDataB.density.mean - mapDataA.density.mean;

maxDelta = max(abs(min(meanDelta)), abs(max(meanDelta)))


%%

clc;
close all;
fig = 1;

orientation = 'YZ';

xLimsData = (mapDataA.positionGrid(1,1)) / normLength;
yLimsData = ([min(mapDataA.positionGrid(:,2)); max(mapDataA.positionGrid(:,2))]) / normLength;
zLimsData = ([min(mapDataA.positionGrid(:,3)); max(mapDataA.positionGrid(:,3))]) / normLength;

xLimsPlot = [0.3; 4.6257662];
yLimsPlot = [-0.5; 0.5];
zLimsPlot = [0; 0.5];

spatialRes = 0.5e-3;
    
positionData = mapDataA.positionGrid;
scalarData = meanDelta;
mapPerim = [];
nPlanes = 1;
planeNo = 1;
cMap = cool2warm(32);
refPoint = [];
figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes

figName = ['Delta_', caseID_A, '_v_', caseID_B];
contourlines = [];
figSubtitle = ' ';
cLims = [-0.16; 0.16]; % Inter-Plane Comparisons at +-0.45
                       % 1.0L Comparisons at +-0.24
                       % 1.5L Comparisons at +-0.32
                       % 2.0L Comparisons at +-0.16

[fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                       xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                       planeNo, fig, figName, cMap, geometry, contourlines, ...
                                       refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                       zLimsPlot, true, figSave);