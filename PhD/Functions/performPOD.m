%% Snapshot POD Calculator v1.0
% ----
% Performs Shapshot POD on Fluctuating Scalar or Vector Fields
%
% Weiss, Julien: A Tutorial on the Proper Orthogonal Decomposition. In: 2019 AIAA Aviation Forum. 17–21
% June 2019, Dallas, Texas, United States.
% ----
% Usage: PODdata = performPOD()


%% Changelog

% v1.0 - Initial Commit


%% Supported Field Types

% Fluctuating Scalar Field: 'scalar'
% Fluctuating Vector Field: 'vector'


%% Main Function

function [fig, PODdata, Ns, Nt, modesEnergetic, modes80percent] = performPOD(fig, PODdata, PODvar, fieldType, location)

    disp('    Performing POD Using the Snapshot Method...');
    
    Ns = height(PODdata.positionGrid); % Number of Spatial Points
    Nt = height(PODdata.time); % Number of Time Instances
    
    % Initialise Progress Bar
    wB = waitbar(0, 'Assembling Snapshot Matrix', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    
    % Assemble Snapshot Matrix
    if strcmp(fieldType, 'scalar')
        PODdata.snapshotMatrix = zeros(Nt,Ns);
        
        for i = 1:Nt
            PODdata.snapshotMatrix(i,:) = PODdata.(PODvar).prime{i};
            
            waitbar((i / Nt), wB);
        end

    elseif strcmp(fieldType, 'vector')
        uSnapshotMatrix = zeros(Nt,Ns);
        vSnapshotMatrix = uSnapshotMatrix;
        wSnapshotMatrix = uSnapshotMatrix;
        
        for i = 1:Nt
            uSnapshotMatrix(i,:) = PODdata.(PODvar{1}).prime{i};
            vSnapshotMatrix(i,:) = PODdata.(PODvar{2}).prime{i};
            wSnapshotMatrix(i,:) = PODdata.(PODvar{3}).prime{i};
            
            waitbar((i / Nt), wB);
        end

        PODdata.snapshotMatrix = [uSnapshotMatrix, vSnapshotMatrix, wSnapshotMatrix];
    else
        error('Unexpected Field Type');
    end
    
    delete(wB);
    
    % Generate Correlation Matrix
    PODdata.C = (PODdata.snapshotMatrix * PODdata.snapshotMatrix') / (Nt - 1);
    
    % Solve Eigenvalue Problem
    [PODdata.A_mode, PODdata.lambda] = eig(PODdata.C, 'vector');
    
    % Sort Eigenvalues and Eigenvalues in Descending Order
    [PODdata.lambda, index] = sort(PODdata.lambda, 'descend');
    PODdata.A_mode = PODdata.A_mode(:,index); % Temporal Modes
    
    % Calculate Spatial Coefficients
    PODdata.phi_coeff = PODdata.snapshotMatrix' * PODdata.A_mode;
    
    % Normalisation to Match Direct Method
    PODdata.phi_mode = normc(PODdata.phi_coeff); % Spatial Modes
    PODdata.A_coeff = PODdata.snapshotMatrix * PODdata.phi_mode; % Temporal Coefficients
    
    % Identify Mode Energy Content
    PODdata.modeEnergy = (PODdata.lambda / sum(PODdata.lambda)) * 100;
    modesEnergetic = height(find(PODdata.modeEnergy > 1));
    modes80percent = find(cumsum(PODdata.modeEnergy) > 80, 1);
    
    disp(' ');
    
    disp(['    First ', num2str(modesEnergetic), ' Modes Each Contain Greater Than 1% of Total Energy']);
    disp(['    First ', num2str(modes80percent), ' Modes Contain Approximately 80% of Total Energy']);
    
    % Figure Setup
    fig = fig + 1;
    
    figName = [location, '_', PODvar, '_POD_Energy_Content'];
    
    set(figure(fig), 'outerPosition', [25, 25, 1275, 850], 'name', figName);
    set(gca, 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
             'fontSize', 20, 'layer', 'top');
    hold on;
    
    % Plot
    plot(PODdata.modeEnergy(1:((ceil(modesEnergetic / 10) * 10) - 1)), 'lineWidth', 1.5, 'marker', 'o', 'color', ([74, 24, 99] / 255));
    
    % Figure Formatting
    axis on;
    box on;
    grid off;
    xlim([0; (ceil(modesEnergetic / 10) * 10)]);
    ylim([0; (ceil(max(PODdata.modeEnergy)/10) * 10)]);
    tickData = (0:(((ceil(modesEnergetic / 10) * 10) - 0) / 5):(ceil(modesEnergetic / 10) * 10));
    xticks(tickData(2:(end - 1)));
    tickData = (0:(((ceil(max(PODdata.modeEnergy)/10) * 10) - 0) / 5):(ceil(max(PODdata.modeEnergy)/10) * 10));
    yticks(tickData(2:(end - 1)));
    xlabel({' ', '{\bf{Mode}}'}, 'fontName', 'LM Roman 12');
    ylabel({'{\bf{Energy Content (\it{%})}}', ' '}, 'fontName', 'LM Roman 12');
    set(gca, 'outerPosition', [0.05, 0.05, 0.9, 0.9]);
    hold off;
    
    pause(2);
    exportgraphics(gcf, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);

end