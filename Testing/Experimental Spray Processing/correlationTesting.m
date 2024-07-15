run preamble;

figSave = false; % Save .fig File(s);

nModes = 1800;

%%

% PODdataA = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayPOD/Far_Field_Soiling_07_22/SB_1.0L_120s_15Hz_02/density/T0067_T120000_F15.mat', 'PODdata').PODdata;
% PODdataB = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayPOD/Far_Field_Soiling_07_22/SB_1.0L_120s_15Hz_01/density/T0067_T120000_F15.mat', 'PODdata').PODdata;

% PODdataA = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayPOD/Far_Field_Soiling_07_22/ST_1.0L_120s_15Hz_01/density/T0067_T120000_F15.mat', 'PODdata').PODdata;
% PODdataB = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayPOD/Far_Field_Soiling_07_22/ST_1.0L_120s_15Hz_02/density/T0067_T120000_F15.mat', 'PODdata').PODdata;

PODdataA = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayPOD/Far_Field_Soiling_07_22/RSST_1.0L_120s_15Hz_02/density/T0067_T120000_F15.mat', 'PODdata').PODdata;
PODdataB = load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayPOD/Far_Field_Soiling_07_22/RSST_1.0L_120s_15Hz_01/density/T0067_T120000_F15.mat', 'PODdata').PODdata;


%%

r = zeros([nModes,1]);
rMode = r;

for i = 1:nModes
    phiA = rescale(PODdataA.POD.phi(:,i), -1, 1);
    
    for j = (i - 5):(i + 5)
        
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

r_AB = [rMode, r]; clear rMode r;

disp(' ');
disp(' ');


%%

p = polyfit(log(1:nModes), log(r_AB(:,2)), 1);
rFit = exp(polyval(p, log(1:0.25:nModes)));

% Initialise Figure
fig = fig + 1;
figName = 'Mode_Correlation_AB';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot Modes
scatter((1:nModes), r_AB(:,2), 10, 'markerFaceColor', graphColours(1), 'markerEdgeColor', graphColours(1));
plot((1:0.25:nModes), rFit, 'color', graphColours(2), 'lineWidth', 2)

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([0; 100]);
ylim([0; 1]);
tickData = (20:20:80);
xticks(tickData);
tickData = (0.2:0.2:0.8);
yticks(tickData);
xlabel({'{Mode}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{Maximum Correlation}'}, 'interpreter', 'latex');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');