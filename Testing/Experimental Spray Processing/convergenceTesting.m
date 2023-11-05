clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

%%%

% load('/mnt/Processing/Data/Experimental/MATLAB/planarContaminantMap/SB_1.0L_120s_15Hz_01/T0067_T120000_F15_Norm.mat');
load('/mnt/Processing/Data/Experimental/MATLAB/planarContaminantMap/RSST_1.0L_600s_03Hz_01/T0333_T600000_F3_Norm.mat');
nTimes = height(mapData.time);

%%

tic;

%%%

close all; 

growingMean = mapData.density.inst;

for i = 1:nTimes
    
    if i == 1
        continue;
    else
        growingMean{i} = growingMean{i-1} + mapData.density.inst{i};
    end
    
end

for i = 1:nTimes
    growingMean{i} = growingMean{i} / i;
end

%%%

growingMeanCorr = zeros([height(nTimes),1]);

for i = 1:nTimes
    growingMeanCorr(i) = corr(growingMean{i}, mapData.density.mean);
end

%%%

growingMeanMean = zeros([height(nTimes),1]);

for i = 1:nTimes
    growingMeanMean(i) = mean(growingMean{i});
end

%%%

toc;

%%

figure();
hold on;
plot(growingMeanMean);
yline(0.99 * growingMeanMean(nTimes));
yline(1.01  * growingMeanMean(nTimes));
set(gca, 'xScale', 'log');
hold off;

% figure();
% hold on;
% plot(growingMeanCorr);
% set(gca, 'xScale', 'log');
% hold off;

% contourlines = [];
% cLims = [0; max(cellfun(@max, growingMeanTest))];
% 
% figHold = fig;
% 
% for i = 100:100:nTimes
% 
%     if i ~= 100
%         clf(fig);
%         fig = figHold;
%     end
% 
%     scalarData = growingMeanTest{i};
%     figTime = num2str(mapData.time(i), '%.3f');
%     figName = ['Growing_Mean_T', erase(figTime, '.')];
%     figSubtitle = [figTime, ' \it{s}'];
% 
%     [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
%                                            xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
%                                            planeNo, fig, figName, cMap, geometry, contourlines, ...
%                                            refPoint, figTitle, figSubtitle, cLims, xLimsPlot, ...
%                                            yLimsPlot, zLimsPlot, normDims, figSave);
% end