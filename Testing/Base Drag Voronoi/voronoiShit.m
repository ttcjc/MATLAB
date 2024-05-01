CLA;

cd ~/Downloads/'Base Drag Voronoi'/;

% % load Exp_SB.mat;
% load CFD_SB.mat;
% % load CFD_SB_uncorrected.mat;
% 
% BoundingBox = [
%                -0.1945, +0.0500
%                -0.1945, +0.3390
%                +0.1945, +0.3390
%                +0.1945, +0.0500
%                -0.1945, +0.0500
%               ];

% % load Exp_ST.mat;
% load CFD_ST.mat;
% 
% BoundingBox = [
%                -0.1945 + (0.045 * sind(20)), +0.0500
%                -0.1945 + (0.045 * sind(20)), +0.3390
%                +0.1945 - (0.045 * sind(20)), +0.3390
%                +0.1945 - (0.045 * sind(20)), +0.0500
%                -0.1945 + (0.045 * sind(20)), +0.0500
%               ];

% % load Exp_RSST.mat;
% load CFD_RSST.mat;
% 
% BoundingBox = [-0.1945,                      +0.0500
%                -0.1945,                      +0.1945
%                -0.1945 + (0.045 * sind(16)), +0.1945
%                -0.1945 + (0.045 * sind(16)), +0.3390
%                +0.1945 - (0.045 * sind(16)), +0.3390
%                +0.1945 - (0.045 * sind(16)), +0.1945
%                +0.1945,                      +0.1945
%                +0.1945,                      +0.0500
%                -0.1945,                      +0.0500
%               ];

% load CFD_FS.mat;
% load CFD_FS_uncorrected.mat;
load CFD_FS_coupled.mat;
% load CFD_FS_halfTread.mat;
% load CFD_FS_20deg.mat;

BoundingBox = [
               -0.778, +0.1820
               -0.778, +1.3380
               +0.778, +1.3380
               +0.778, +0.1820
               -0.778, +0.1820
              ];


TappingCoords = pData.positionGrid(:,[2,3]);

PressureData = pData.Cp.mean; PressureData(isnan(PressureData)) = 0;

FigNo = 1;

Range = 'auto';

FrontOrBase = true;

PlotReq = false;

plotType = 'main';

tileSelect = [];

tic;
[AverageCPVoronoi, FigNo] = VoronoiCDCalcV3(TappingCoords, BoundingBox, PressureData, FigNo, Range, ...
                                            FrontOrBase, PlotReq, plotType, tileSelect);
executionTime = toc;

baseCd = round(abs(AverageCPVoronoi), 3);

disp(['Voronoi Estimate of Base Cd: ', num2str(baseCd)]);
disp(['    Run Time: ', num2str(toc), ' s']);
                                        