%% Snapshot POD Calculator v1.2
% ----
% Performs Shapshot POD on Fluctuating Scalar or Vector Fields
%
% J. Weiss
% "A Tutorial on the Proper Orthogonal Decomposition"
% 2019 AIAA Aviation Forum, 17-21 June 2019, Dallas, Texas, United States
% ----
% Usage: [fig, PODdata, modesEnergetic, modes80percent, Ns, Nt] = performPOD(fig, PODdata, field, ...
%                                                                            fieldType, location, figSave)
%
%        'fig'       -> Figure Number
%        'PODdata'   -> Structure Containing Position and Field Data
%        'field'     -> Field Variable Used to Perform POD Stored as String
%        'fieldType' -> Desired Field Type Stored as String
%        'location'  -> Data Location Identifier Stored as String
%        'figSave'     -> Save .fig File [True/False]


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Formatting Updates
% v1.2 - Further Formatting Updates

%% Supported Field Types

% Fluctuating Scalar Field: 'scalar'
% Fluctuating Vector Field: 'vector'


%% Main Function

function [fig, PODdata, modesEnergetic, modes90, Ns, Nt] = performPOD(fig, PODdata, field, ...
                                                                             fieldType, location, figSave)

    PODdata.POD.field = field;
    
    Ns = height(PODdata.positionGrid); % Number of Spatial Points
    Nt = height(PODdata.time); % Number of Time Instances

    disp('    Assembling Snapshot Matrix...');
    
    % Initialise Progress Bar
    wB = waitbar(0, 'Assembling Snapshot Matrix', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    
    % Assemble Snapshot Matrix
    switch fieldType

        case 'scalar'
            snapshotMatrix = zeros(Nt,Ns);
            
            for i = 1:Nt
                snapshotMatrix(i,:) = PODdata.(field).prime{i};
                
                % Update Waitbar
                waitbar((i / Nt), wB);
            end
            clear i;

        case 'vector'
            uSnapshotMatrix = zeros(Nt,Ns);
            vSnapshotMatrix = uSnapshotMatrix;
            wSnapshotMatrix = uSnapshotMatrix;
            
            for i = 1:Nt
                uSnapshotMatrix(i,:) = PODdata.(field{1}).prime{i};
                vSnapshotMatrix(i,:) = PODdata.(field{2}).prime{i};
                wSnapshotMatrix(i,:) = PODdata.(field{3}).prime{i};
                
                % Update Waitbar
                waitbar((i / Nt), wB);
            end
            clear i;
    
            snapshotMatrix = [uSnapshotMatrix, vSnapshotMatrix, wSnapshotMatrix];
    
    end
    
    delete(wB);

    disp(' ');

    disp('    Performing POD Using the Snapshot Method...');
    
    % Generate Correlation Matrix
    C = (snapshotMatrix * snapshotMatrix') / (Nt - 1);
    
    % Solve Eigenvalue Problem
    [alpha_s, lambda] = eig(C, 'vector');
    
    % Sort Eigenvalues and Eigenvalues in Descending Order
    [lambda, index] = sort(lambda, 'descend');
    alpha_s = alpha_s(:,index); % Temporal Modes
    
    % Calculate Spatial Coefficients
    phi_s = snapshotMatrix' * alpha_s;
    
    % Normalisation to Match Direct Method
    PODdata.POD.phi = normc(phi_s); % Spatial Modes
    PODdata.POD.alpha = snapshotMatrix * PODdata.POD.phi; % Temporal Coefficients
    
    % Identify Mode Energy Content
    PODdata.POD.modeEnergy = (lambda / sum(lambda)) * 100;
    modesEnergetic = height(find(PODdata.POD.modeEnergy > 1));
    modes90 = find(cumsum(PODdata.POD.modeEnergy) > 90, 1);
    
    disp(['        The First ', num2str(modesEnergetic), ' Modes Each Contain Greater Than 1% of the Total Captured Energy Content']);
    disp(['        The First ', num2str(modes90), ' Modes Contain Approximately 90% of the Total Captured Energy Content']);
    
    % Initialise Figure
    fig = fig + 1;
    
    switch fieldType

        case 'scalar'
            figName = [location, '_POD_', field, '_Mode_Energy_Content'];

        case 'vector'
            figName = [location, '_POD_', cell2mat(field), '_Mode_Energy_Content'];

    end
    
    set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                     'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
    pause(0.5);
    hold on;
    set(gca, 'positionConstraint', 'outerPosition', ...
             'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 20, 'layer', 'top');
    
    % Plot Mode Energies
    plot(PODdata.POD.modeEnergy(1:49), 'color', graphColours(1), 'lineWidth', 2)
    scatter((1:49), PODdata.POD.modeEnergy(1:49), 30, 'markerFaceColor', graphColours(1), 'markerEdgeColor', graphColours(1));
    
    % Format Figure
    title('{-----}', 'interpreter', 'latex');
    subtitle('{ }');
    axis on;
    box on;
    grid off;
    xlim([0; 50]);
    ylim([0; 16]);
    tickData = (10:10:40);
    xticks(tickData);
    tickData = (3:3:12);
    yticks(tickData);
    xlabel({'{Mode}'; '{-----}'}, 'interpreter', 'latex');
    ylabel({'{-----}'; '{Energy Content (\%)}'}, 'interpreter', 'latex');
    tightInset = get(gca, 'TightInset');
    set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                               (tightInset(2) + 0.00625), ...
                               (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                               (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
    pause(0.5);
    hold off;

    % Save Figure
    print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');

    if figSave
        savefig(gcf, [userpath, '/Output/Figures/', figName, '.fig']);
    end
    
end