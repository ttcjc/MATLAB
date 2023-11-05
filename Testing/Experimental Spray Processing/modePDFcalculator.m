% clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

figSave = false; % Save .fig File(s);


%%

% PODdata = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayPOD/SB_1.0L_120s_15Hz_02/density/T0067_T120000_F15_Norm.mat', 'PODdata').PODdata;
% PODdata = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayPOD/ST_1.0L_120s_15Hz_01/density/T0067_T120000_F15_Norm.mat', 'PODdata').PODdata;
% PODdata = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayPOD/RSST_1.0L_120s_15Hz_02/density/T0067_T120000_F15_Norm.mat', 'PODdata').PODdata;

% PODdata = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayPOD/SB_1.0L_600s_03Hz_01/density/T0333_T600000_F3_Norm.mat', 'PODdata').PODdata;


%%

% modeOfInterest = 1;
% 
% alphaScaled = rescale(PODdata.POD.alpha(:,modeOfInterest), -1, 1);
% 
% alphaFit = fitdist(alphaScaled, 'kernel');
% alphaPDF = pdf(alphaFit, (-1.2:1e-3:1.2));
% alphaPDF = rescale(alphaPDF, 0, 1);
% 
% % Initialise Figure
% fig = fig + 1;
% figName = ['Alpha_PDF_', caseID];
% set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
%                  'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
% pause(0.5);
% hold on;
% set(gca, 'positionConstraint', 'outerPosition', ...
%          'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 20, 'layer', 'top');
% 
% % Plot
% plot((-1.2:1e-3:1.2), alphaPDF, 'color', graphColours(1), 'lineWidth', 2);
% plot([0; 0], [0; 1.2], 'color', [0.15 0.15 0.15], 'lineWidth', 2);
% 
% % Format Figure
% title('{-----}', 'interpreter', 'latex');
% subtitle('{ }');
% axis on;
% box on;
% grid off;
% xlim([-1.2; 1.2]);
% ylim([0; 1.2]);
% tickData = (-0.72:0.48:0.72);
% xticks(tickData);
% tickData = [];
% yticks(tickData);
% xtickformat('%+.2g');
% xlabel({'{$\alpha_{n}$}'; '{-----}'}, 'interpreter', 'latex');
% ylabel({'{-----}'; '{Probability Density}'}, 'interpreter', 'latex');
% tightInset = get(gca, 'TightInset');
% set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
%                            (tightInset(2) + 0.00625), ...
%                            (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
%                            (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
% pause(0.5);
% hold off;
% 
% print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');


%%

PDFdata.PODdataA = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayPOD/SB_1.0L_120s_15Hz_03/density/T0067_T120000_F15_Norm.mat', 'PODdata').PODdata;
PDFdata.PODdataB = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayPOD/ST_1.0L_120s_15Hz_03/density/T0067_T120000_F15_Norm.mat', 'PODdata').PODdata;
PDFdata.PODdataC = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayPOD/RSST_1.0L_120s_15Hz_03/density/T0067_T120000_F15_Norm.mat', 'PODdata').PODdata;
% 
setsAll = fieldnames(PDFdata);

%%

modeOfInterest = 3;

% Initialise Figure
fig = fig + 1;
figName = ['Alpha_', num2str(modeOfInterest), '_PDF'];
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 20, 'layer', 'top');

% Calculate & Plot PDF
for i = 1:height(setsAll)
    setsCurrent = setsAll{i};
    
    alphaScaled = rescale(PDFdata.(setsCurrent).POD.alpha(:,modeOfInterest), -1, 1);

    alphaFit = fitdist(alphaScaled, 'kernel');
    alphaPDF = pdf(alphaFit, (-1.2:1e-3:1.2));
    alphaPDF = rescale(alphaPDF, 0, 1);
    
    plot((-1.2:1e-3:1.2), alphaPDF, 'color', graphColours(i), 'lineWidth', 2);
end
clear i;

% Plot Zero Reference
plot([0; 0], [0; 1.2], 'color', [0.15 0.15 0.15], 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([-1.2; 1.2]);
ylim([0; 1.2]);
tickData = (-0.72:0.48:0.72);
xticks(tickData);
tickData = [];
yticks(tickData);
xtickformat('%+.2g');
xlabel({'{$\alpha_{n}$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{Probability Density}'}, 'interpreter', 'latex');
legend({'\textit{Config A}', '\textit{Config B}', '\textit{Config C}'}, ...
       'location', 'northEast', 'fontSize', 16, 'interpreter', 'latex'); legend('boxoff');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');