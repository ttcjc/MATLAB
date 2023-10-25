%% Snapshot POD Calculator v1.2
% ----
% Performs Shapshot POD on Fluctuating Scalar or Vector Fields
%
% J. Weiss
% "A Tutorial on the Proper Orthogonal Decomposition"
% 2019 AIAA Aviation Forum, 17-21 June 2019, Dallas, Texas, United States
% ----
% Usage: [fig, PODdata, modesEnergetic, modes80percent, Ns, Nt] = performPOD(fig, PODdata, PODvar, ...
%                                                                            fieldType, location, figSave)
%
%        'fig'       -> Figure Number
%        'PODdata'   -> Structure Containing Position and Field Data
%        'PODvar'    -> Field Variable Used to Perform POD Stored as String
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

function [fig, PODdata, modesEnergetic, modes80percent, Ns, Nt] = performPOD(fig, PODdata, PODvar, ...
                                                                             fieldType, location, figSave)

    Ns = height(PODdata.positionGrid); % Number of Spatial Points
    Nt = height(PODdata.time); % Number of Time Instances

    disp('    Assembling Snapshot Matrix...');
    
    % Initialise Progress Bar
    wB = waitbar(0, 'Assembling Snapshot Matrix', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    
    % Assemble Snapshot Matrix
    switch fieldType

        case 'scalar'
            PODdata.snapshotMatrix = zeros(Nt,Ns);
            
            for i = 1:Nt
                PODdata.snapshotMatrix(i,:) = PODdata.(PODvar).prime{i};
                
                % Update Waitbar
                waitbar((i / Nt), wB);
            end
            clear i;

        case 'vector'
            uSnapshotMatrix = zeros(Nt,Ns);
            vSnapshotMatrix = uSnapshotMatrix;
            wSnapshotMatrix = uSnapshotMatrix;
            
            for i = 1:Nt
                uSnapshotMatrix(i,:) = PODdata.(PODvar{1}).prime{i};
                vSnapshotMatrix(i,:) = PODdata.(PODvar{2}).prime{i};
                wSnapshotMatrix(i,:) = PODdata.(PODvar{3}).prime{i};
                
                % Update Waitbar
                waitbar((i / Nt), wB);
            end
            clear i;
    
            PODdata.snapshotMatrix = [uSnapshotMatrix, vSnapshotMatrix, wSnapshotMatrix];
    
    end
    
    delete(wB);

    disp(' ');

    disp('    Performing POD Using the Snapshot Method...');
    
    % Generate Correlation Matrix
    PODdata.C = (PODdata.snapshotMatrix * PODdata.snapshotMatrix') / (Nt - 1);
    
    % Solve Eigenvalue Problem
    [A_mode, lambda] = eig(PODdata.C, 'vector');
    
    % Sort Eigenvalues and Eigenvalues in Descending Order
    [lambda, index] = sort(lambda, 'descend');
    A_mode = A_mode(:,index); % Temporal Modes
    
    % Calculate Spatial Coefficients
    phi_coeff = PODdata.snapshotMatrix' * A_mode;
    
    % Normalisation to Match Direct Method
    PODdata.phi_mode = normc(phi_coeff); % Spatial Modes
    PODdata.A_coeff = PODdata.snapshotMatrix * PODdata.phi_mode; % Temporal Coefficients
    
    % Identify Mode Energy Content
    PODdata.modeEnergy = (lambda / sum(lambda)) * 100;
    modesEnergetic = height(find(PODdata.modeEnergy > 1));
    modes80percent = find(cumsum(PODdata.modeEnergy) > 80, 1);
    
    disp(' ');
    
    disp(['    First ', num2str(modesEnergetic), ' Modes Each Contain Greater Than 1% of Total Energy']);
    disp(['    First ', num2str(modes80percent), ' Modes Contain Approximately 90% of Total Energy']);
    
    % Initialise Figure
    fig = fig + 1;
    
    switch fieldType

        case 'scalar'
            figName = [location, '_POD_', PODvar, '_Mode_Energy_Content'];

        case 'vector'
            figName = [location, '_POD_', cell2mat(PODvar), '_Mode_Energy_Content'];

    end
    
    set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                     'units', 'pixels', 'outerPosition', [50, 50, 795, 880])
    set(gca, 'positionConstraint', 'outerPosition', ...
             'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 20, 'layer', 'top');
    hold on;
    
    % Plot Modes
    bar(PODdata.modeEnergy, 0.75, ...
        'lineWidth', 2, 'faceColor', graphColours(1));
    
    figTitle = '.';
    figSubtitle = ' ';
    
    % Format Figure
    title(figTitle, 'color', ([254, 254, 254] / 255));
    subtitle(figSubtitle);
    axis on;
    box on;
    grid off;
    xlim([0; 50]);
    ylim([0; 20]);
    tickData = (10:10:40);
    xticks(tickData);
    tickData = (4:4:16);
    yticks(tickData);
    xlabel('{Mode}', 'interpreter', 'latex');
    ylabel('{Energy Content (\%)}', 'interpreter', 'latex');
    hold off;

    tightInset = get(gca, 'TightInset');
    set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                               (tightInset(2) + 0.00625), ...
                               (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                               (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
    hold off;

    pause(1);

    % Save Figure
    print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');

    if figSave
        savefig(gcf, [userpath, '/Output/Figures/', figName, '.fig']);
    end
    
end
