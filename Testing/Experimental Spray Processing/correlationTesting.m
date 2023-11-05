clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

figSave = false; % Save .fig File(s);

nModes = 1800;

%%

PODdataA = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayPOD/SB_1.0L_120s_15Hz_01/density/T0067_T120000_F15_Norm.mat', 'PODdata').PODdata;
PODdataB = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayPOD/SB_1.0L_120s_15Hz_02/density/T0067_T120000_F15_Norm.mat', 'PODdata').PODdata;
PODdataC = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayPOD/SB_1.0L_120s_15Hz_03/density/T0067_T120000_F15_Norm.mat', 'PODdata').PODdata;


%%

r = zeros([nModes,1]);
rMode = r;

tic;
for i = 1:nModes
    phiA = rescale(PODdataA.POD.phi(:,i), -1, 1);
    
    for j = (i - 10):(i + 10)
        
        if j < 1 || j > 1800
            continue;
        else
            phiB = rescale(PODdataB.POD.phi(:,j), -1, 1);
            rTemp = abs(corr(phiA, phiB));
        
            if max(abs(rTemp), abs(r(i))) == abs(rTemp)
                r(i) = rTemp;
                rMode(i) = j;
            end
            
        end
        
    end
    clear j rTemp;
    
end
clear i;
toc;

r_AB = [rMode, r]; clear rMode r;

disp(' ');
disp(' ');


%%

r = zeros([nModes,1]);
rMode = r;

tic;
for i = 1:nModes
    phiA = rescale(PODdataA.POD.phi(:,i), -1, 1);
    
    for j = (i - 10):(i + 10)
        
        if j < 1 || j > 1800
            continue;
        else
            phiB = rescale(PODdataC.POD.phi(:,j), -1, 1);
            rTemp = abs(corr(phiA, phiB));
        
            if max(abs(rTemp), abs(r(i))) == abs(rTemp)
                r(i) = rTemp;
                rMode(i) = j;
            end
            
        end
        
    end
    clear j rTemp;
    
end
clear i;
toc;

r_AC = [rMode, r]; clear rMode r;


%%

p = polyfit(log(1:nModes), log(r_AB(:,2)), 1);
rFit = exp(polyval(p, log(1:0.25:nModes)));

% Initialise Figure
fig = fig + 1;
set(figure(fig), 'name', 'Mode_Correlation_AB', 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880])
set(gca, 'positionConstraint', 'outerPosition', ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 20, 'layer', 'top');
hold on;

% Plot Modes
scatter((1:nModes), r_AB(:,2), 10, 'markerFaceColor', graphColours(1), 'markerEdgeColor', graphColours(1));
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


p = polyfit(log(1:nModes), log(r_AC(:,2)), 1);
rFit = exp(polyval(p, log(1:0.25:nModes)));

% Initialise Figure
fig = fig + 1;
set(figure(fig), 'name', 'Mode_Correlation_AC', 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880])
set(gca, 'positionConstraint', 'outerPosition', ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 20, 'layer', 'top');
hold on;

% Plot Modes
scatter((1:nModes), r_AC(:,2), 10, 'markerFaceColor', graphColours(1), 'markerEdgeColor', graphColours(1));
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