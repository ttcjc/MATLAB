run preamble;

figSave = false; % Save .fig File(s);


%%

PODdataA = load('/mnt/Processing/Data/Experimental/MATLAB/planarContaminantPOD/SB_1.0L_120s_15Hz_03/density/T0067_T120000_F15_Norm.mat', 'PODdata').PODdata;
PODdataB = load('/mnt/Processing/Data/Experimental/MATLAB/planarContaminantPOD/SB_1.0L_120s_15Hz_02/density/T0067_T120000_F15_Norm.mat', 'PODdata').PODdata;

%%

clc;

nModes = 1800;

r = zeros([nModes,1]);
rMode = r;

for i = 1:nModes
    phiA = rescale(PODdataA.phi_mode(:,i), -1, 1);
    
    for j = 1:nModes
        phiB = rescale(PODdataB.phi_mode(:,j), -1, 1);
        rTemp = abs(corr(phiA, phiB));
        
        if max(abs(rTemp), abs(r(i))) == abs(rTemp)
            r(i) = rTemp;
            rMode(i) = j;
        end
        
    end
    clear rTemp;
    
end

%%%

p = polyfit(log(1:nModes), log(r), 1);
rFit = exp(polyval(p, log(1:0.25:nModes)));

% Initialise Figure
fig = fig + 1;
set(figure(fig), 'name', 'POD Filtering', 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880])
set(gca, 'positionConstraint', 'outerPosition', ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 20, 'layer', 'top');
hold on;

% Plot Modes
scatter((1:nModes), r, 10, 'markerFaceColor', graphColours(1), 'markerEdgeColor', graphColours(1));
plot((1:0.25:nModes), rFit, 'color', graphColours(2), 'lineWidth', 2)

figTitle = '.';
figSubtitle = ' ';

% Format Figure
title(figTitle, 'color', ([254, 254, 254] / 255));
subtitle(figSubtitle);
axis on;
box on;
grid off;
xlim([1; nModes]);
ylim([0; 1]);
tickData = (300:300:1500);
xticks(tickData);
tickData = (0.2:0.2:0.8);
yticks(tickData);
xlabel('{Mode}', 'interpreter', 'latex');
ylabel('{Maximum Correlation}', 'interpreter', 'latex');
hold off;

tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
hold off;

pause(1);