run preamble;

figSave = false; % Save .fig File(s);


%%

PDFdata.PODdataA = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayPOD/Far_Field_Soiling_07_22/SB_1.0L_120s_15Hz_02/density/T0067_T120000_F15.mat', 'PODdata').PODdata;
PDFdata.PODdataB = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayPOD/Far_Field_Soiling_07_22/ST_1.0L_120s_15Hz_01/density/T0067_T120000_F15.mat', 'PODdata').PODdata;
PDFdata.PODdataC = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayPOD/Far_Field_Soiling_07_22/RSST_1.0L_120s_15Hz_02/density/T0067_T120000_F15.mat', 'PODdata').PODdata;
% 
setsAll = fieldnames(PDFdata);

%%

modeOfInterest = 1;

% Initialise Figure
fig = fig + 1;
figName = ['Alpha_', num2str(modeOfInterest), '_PDF'];
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Calculate & Plot PDF
maxVal = 0;

for i = 1:height(setsAll)
    setsCurrent = setsAll{i};
    
    alphaScaled = PDFdata.(setsCurrent).POD.alpha(:,modeOfInterest) / max(abs(PDFdata.(setsCurrent).POD.alpha(:,modeOfInterest)));

    alphaFit = fitdist(alphaScaled, 'kernel');
    alphaPDF = pdf(alphaFit, (-1.2:1e-3:1.2));
    
    maxVal = max(maxVal, max(alphaPDF));
    
    plot((-1.2:1e-3:1.2), alphaPDF, 'color', graphColours(i), 'lineWidth', 2);
end
clear i;

% Plot Zero Reference
plot([0; 0], [-1e3; 1e3], 'color', [0.15 0.15 0.15], 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([-1.2; 1.2]);
ylim([0; (1.05 * maxVal)]);
tickData = (-0.72:0.48:0.72);
xticks(tickData);
tickData = [];
yticks(tickData);
xtickformat('%+.2g');
xlabel({'{$\alpha_{_{n}}$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{Probability Density}'}, 'interpreter', 'latex');
legend({'\textit{Config~A}', '\textit{Config~B}', '\textit{Config~C}'}, ...
       'location', 'northEast', 'fontSize', 16, 'interpreter', 'latex'); legend('boxoff');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');