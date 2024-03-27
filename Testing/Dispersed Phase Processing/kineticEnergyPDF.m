run preamble;

%#ok<*UNRCH>

caseID = '20deg';

load(['~/MATLAB/Testing/Dispersed Phase Processing/LagDataPlaneRawFS_', caseID, '.mat']);

LagData = orderfields(LagData, [4,1,2,3]);

planes = fieldnames(LagData);

energyRange = linspace(0, 400, 1e3)';

%% Concatenate Instantaneous Data

for i = 1:height(planes)
    LagData.(planes{i}).d = double(cell2mat(LagData.(planes{i}).d));
    LagData.(planes{i}).nParticle = double(cell2mat(LagData.(planes{i}).nParticle));
    LagData.(planes{i}).Up = double(cell2mat(LagData.(planes{i}).Up));
    
    LagData.(planes{i}).Up = sqrt(LagData.(planes{i}).Up(:,1).^2 + ...
                                  LagData.(planes{i}).Up(:,2).^2 + ...
                                  LagData.(planes{i}).Up(:,3).^2);
end
clear i;


%% Calculate Kinetic Energy Distribution

for i = 1:height(planes)
    
    % Calculate Massic Kinetic Energy
    if strcmp(caseID, 'halfTread')
        LagData.(planes{i}).Ek = (LagData.(planes{i}).nParticle .* ...
                                 (0.5 * (1000 * ((1 / 12) * tau * LagData.(planes{i}).d.^3)) .* ...
                                 (LagData.(planes{i}).Up.^2))) / 5.884562916496705e-06;
    else
        LagData.(planes{i}).Ek = (LagData.(planes{i}).nParticle .* ...
                                 (0.5 * (1000 * ((1 / 12) * tau * LagData.(planes{i}).d.^3)) .* ...
                                 (LagData.(planes{i}).Up.^2))) / 1.079980344450543e-05;
    end

    % Calculate PDF
    Ek_Fit = fitdist(LagData.(planes{i}).Ek, 'kernel');
    Ek_PDF = pdf(Ek_Fit, energyRange);
    
    % Initialise Figure
    if i == 1
        fig = fig + 1;
        figName = ['Kinetic_Energy_Distribution_', caseID];
        set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                     'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
        pause(0.5);
        hold on;
        set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
                 'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');
    end
    
    % Plot Kinetic Energy Distribution
    plot(energyRange, Ek_PDF, 'color', graphColours(i), 'lineWidth', 2);
    
    % Format Figure
    if i == height(planes)
        title('{-----}', 'interpreter', 'latex');
        subtitle('{ }');
        axis on;
        box on;
        grid off;
        xlim([0; 400]);
        ylim([0; 0.025]);
        tickData = (80:80:320);
        xticks(tickData);
        tickData = (5e-3:5e-3:20e-3);
        yticks(tickData);
        currentAxis = gca; currentAxis.YAxis.Exponent = -3;
        xlabel({'{Massic Kinetic Energy ($m^{2} {\cdot} s^{-2}$)}'; '{-----}'}, 'interpreter', 'latex');
        ylabel({'{-----}'; '{Probability Density}'}, 'interpreter', 'latex');
        legend({'$1\,\ell$', ...
                '$2\,\ell$', ...
                '$3\,\ell$', ...
                '$4\,\ell$'}, 'location', 'northEast', 'orientation', 'vertical', 'interpreter', 'latex', ...
               'fontSize', 18, 'box', 'off');
        tightInset = get(gca, 'TightInset');
        set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                                   (tightInset(2) + 0.00625), ...
                                   (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                                   (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
        pause(0.5);
        hold off;

        % Save Figure
        print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');
    end
    
end
clear i;