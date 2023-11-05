clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

figSave = false; % Save .fig File(s);

%%%

% % Select Relevant Geometry and Define Bounding Box
% [geometry, xDims, yDims, zDims, spacePrecision, normDims, normLength] = selectGeometry(true);


%%

clc;
close all;


%%

% 1.0L Config Repeatability
% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/SB_1.0L_120s_15Hz_02/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/SB_1.0L_120s_15Hz_01/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/SB_1.0L_120s_15Hz_03/T0067_T120000_F15_Norm.mat';

% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/ST_1.0L_120s_15Hz_01/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/ST_1.0L_120s_15Hz_02/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/ST_1.0L_120s_15Hz_03/T0067_T120000_F15_Norm.mat';

% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/RSST_1.0L_120s_15Hz_02/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/RSST_1.0L_120s_15Hz_01/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/RSST_1.0L_120s_15Hz_03/T0067_T120000_F15_Norm.mat';


% Per-Plane Config Variations
% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/SB_1.0L_120s_15Hz_02/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/ST_1.0L_120s_15Hz_01/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/RSST_1.0L_120s_15Hz_02/T0067_T120000_F15_Norm.mat';

% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/SB_1.5L_120s_15Hz_01/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/ST_1.5L_120s_15Hz_02/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/RSST_1.5L_120s_15Hz_03/T0067_T120000_F15_Norm.mat';

% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/SB_2.0L_120s_15Hz_03/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/ST_2.0L_120s_15Hz_03/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/RSST_2.0L_120s_15Hz_02/T0067_T120000_F15_Norm.mat';


% Other Stuff
% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/SB_1.0L_120s_15Hz_02/T0067_T120000_F15_Norm.mat';
% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/SB_1.5L_120s_15Hz_01/T0067_T120000_F15_Norm.mat';
% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/SB_2.0L_120s_15Hz_03/T0067_T120000_F15_Norm.mat';

% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/SB_1.0L_120s_15Hz_02/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/SB_1.5L_120s_15Hz_01/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/SB_2.0L_120s_15Hz_03/T0067_T120000_F15_Norm.mat';


% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/ST_1.0L_120s_15Hz_01/T0067_T120000_F15_Norm.mat';
% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/ST_1.5L_120s_15Hz_02/T0067_T120000_F15_Norm.mat';
% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/ST_2.0L_120s_15Hz_03/T0067_T120000_F15_Norm.mat';

% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/ST_1.0L_120s_15Hz_01/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/ST_1.5L_120s_15Hz_02/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/ST_2.0L_120s_15Hz_03/T0067_T120000_F15_Norm.mat';


% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/RSST_1.0L_120s_15Hz_02/T0067_T120000_F15_Norm.mat';
% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/RSST_1.5L_120s_15Hz_03/T0067_T120000_F15_Norm.mat';
% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/RSST_2.0L_120s_15Hz_02/T0067_T120000_F15_Norm.mat';

% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/RSST_1.0L_120s_15Hz_02/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/RSST_1.5L_120s_15Hz_03/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/RSST_2.0L_120s_15Hz_02/T0067_T120000_F15_Norm.mat';


% Load Data
mapDataA = load(caseA, 'mapData').mapData;
mapDataB = load(caseB, 'mapData').mapData;

caseID_A = load(caseA, 'caseID').caseID;
caseID_B = load(caseB, 'caseID').caseID;


%%

% correlation = corr(mapDataA.density.RMS, mapDataB.density.RMS)

valsA = [min(mapDataA.density.mean), max(mapDataA.density.mean), mean(mapDataA.density.mean(mapDataA.density.mean > 1e-6)), mean(mapDataA.density.RMS(mapDataA.density.mean > 1e-6))]
valsB = [min(mapDataB.density.mean), max(mapDataB.density.mean), mean(mapDataB.density.mean(mapDataB.density.mean > 1e-6)), mean(mapDataB.density.RMS(mapDataB.density.mean > 1e-6))]

% RMSerror = (sqrt((1 / height(mapDataA.time)) * sum((mapDataA.density.mean - mapDataB.density.mean).^2))) / mean(mapDataA.density.mean)


%%z

meanDelta = mapDataB.density.mean - mapDataA.density.mean;

maxDelta = max(abs(min(meanDelta)), abs(max(meanDelta)))


%%

% orientation = 'YZ';
% 
% xLimsData = mapDataA.positionGrid(1,1);
% yLimsData = [min(mapDataA.positionGrid(:,2)); max(mapDataA.positionGrid(:,2))];
% zLimsData = [min(mapDataA.positionGrid(:,3)); max(mapDataA.positionGrid(:,3))];
% 
% xLimsPlot = [0.31875; 2.73575];
% yLimsPlot = [-0.522; 0.522];
% zLimsPlot = [0; 0.522];
% xLimsPlot = round((xLimsPlot / normLength), spacePrecision);
% yLimsPlot = round((yLimsPlot / normLength), spacePrecision);
% zLimsPlot = round((zLimsPlot / normLength), spacePrecision);
% 
% spatialRes = 0.5e-3;
%     
% positionData = mapDataA.positionGrid;
% scalarData = meanDelta;
% mapPerim = [];
% nPlanes = 1;
% planeNo = 1;
% cMap = cool2warm(32);
% refPoint = [];
% figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
% 
% figName = ['Delta_', caseID_A, '_v_', caseID_B];
% contourlines = [];
% figSubtitle = ' ';
% cLims = [-0.16; 0.16]; % Inter-Plane Comparisons at +-0.45
%                        % 1.0L Comparisons at +-0.25
%                        % 1.5L Comparisons at +-0.32
%                        % 2.0L Comparisons at +-0.16
% 
% [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
%                                        xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
%                                        planeNo, fig, figName, cMap, geometry, contourlines, ...
%                                        refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
%                                        zLimsPlot, normDims, figSave);