CLA;

load('/mnt/Processing/Data/Numerical/MATLAB/LagData/Windsor_Upstream_2023/Windsor_SB_wW_Upstream_SC/plane/T12525_T40000_F400_cumulative.mat');
LagData = LagData.X_P1_24625;

load('/mnt/Processing/Data/Numerical/MATLAB/planarVelocityPOD/Windsor_Upstream_2023/Windsor_SB_wW_Upstream_SC/T12525_T40000_F400/X_P1_2462500.mat');


%%

nTimes = height(LagData.time);

massFlux = zeros([nTimes, 1], 'single');

for i = 1:nTimes
    massFlux(i) = sum(LagData.nParticle{i} .* ...
                           (1000 .* ((1 / 12) * tau * LagData.d{i}.^3)));
end
clear i;


%% 

close all;

for i = 1:8
    mass = rescale(massFlux, -1, 1);
    alpha = rescale(PODdata.POD.alpha(:,i), -1, 1);
    
    figure(i);
    set(figure(i), 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
    pause(0.5);
    hold on;
    set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
             'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');
    
    histogram2(mass, alpha, 40, 'displayStyle', 'tile', 'showEmptyBins', 'on');
    
    box on;
    grid off;
    xlim([-1; 1]);
    ylim([-1; 1]);
    xticks(-0.6:0.4:0.6);
    yticks(-0.6:0.4:0.6);
    xtickformat('%+.1f');
    ytickformat('%+.1f');
    xlabel({'{mass}'; '{-----}'}, 'interpreter', 'latex');
    ylabel({'{alpha}'; '{-----}'}, 'interpreter', 'latex');
    tightInset = get(gca, 'TightInset');
    set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                               (tightInset(2) + 0.00625), ...
                               (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                               (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
    pause(0.5);
    hold off;
end