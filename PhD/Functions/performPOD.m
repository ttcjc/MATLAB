%% Snapshot POD Calculator v1.1
% ----
% Performs Shapshot POD on Fluctuating Scalar or Vector Fields
%
% Weiss, Julien: A Tutorial on the Proper Orthogonal Decomposition. In: 2019 AIAA Aviation Forum. 17–21
% June 2019, Dallas, Texas, United States.
% ----
% Usage: [fig, PODdata, modesEnergetic, modes80percent, Ns, Nt] = performPOD(fig, PODdata, PODvar, ...
%                                                                            fieldType, location)
%        'fig'       -> Figure Number
%        'PODdata'   -> Structure Containing Position and Field Data
%        'PODvar'    -> Field Variable Used to Perform POD Stored as String
%        'fieldType' -> Desired Field Type Stored as String
%        'location'  -> Data Location Identifier Stored as String


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Formatting Updates


%% Supported Field Types

% Fluctuating Scalar Field: 'scalar'
% Fluctuating Vector Field: 'vector'


%% Main Function

function [fig, PODdata, modesEnergetic, modes80percent, Ns, Nt] = performPOD(fig, PODdata, PODvar, ...
                                                                             fieldType, location)

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
                snapshotMatrix(i,:) = PODdata.(PODvar).prime{i};
                
                waitbar((i / Nt), wB);
            end

        case 'vector'
            uSnapshotMatrix = zeros(Nt,Ns);
            vSnapshotMatrix = uSnapshotMatrix;
            wSnapshotMatrix = uSnapshotMatrix;
            
            for i = 1:Nt
                uSnapshotMatrix(i,:) = PODdata.(PODvar{1}).prime{i};
                vSnapshotMatrix(i,:) = PODdata.(PODvar{2}).prime{i};
                wSnapshotMatrix(i,:) = PODdata.(PODvar{3}).prime{i};
                
                waitbar((i / Nt), wB);
            end
    
            snapshotMatrix = [uSnapshotMatrix, vSnapshotMatrix, wSnapshotMatrix];
    
    end
    
    delete(wB);

    disp(' ');

    disp('    Performing POD Using the Snapshot Method...');
    
    % Generate Correlation Matrix
    C = (snapshotMatrix * snapshotMatrix') / (Nt - 1);
    
    % Solve Eigenvalue Problem
    [A_mode, lambda] = eig(C, 'vector');
    
    % Sort Eigenvalues and Eigenvalues in Descending Order
    [lambda, index] = sort(lambda, 'descend');
    A_mode = A_mode(:,index); % Temporal Modes
    
    % Calculate Spatial Coefficients
    phi_coeff = snapshotMatrix' * A_mode;
    
    % Normalisation to Match Direct Method
    PODdata.phi_mode = normc(phi_coeff); % Spatial Modes
    PODdata.A_coeff = snapshotMatrix * PODdata.phi_mode; % Temporal Coefficients
    
    % Identify Mode Energy Content
    PODdata.modeEnergy = (lambda / sum(lambda)) * 100;
    modesEnergetic = height(find(PODdata.modeEnergy > 1));
    modes80percent = find(cumsum(PODdata.modeEnergy) > 80, 1);
    
    disp(' ');
    
    disp(['    First ', num2str(modesEnergetic), ' Modes Each Contain Greater Than 1% of Total Energy']);
    disp(['    First ', num2str(modes80percent), ' Modes Contain Approximately 80% of Total Energy']);
    
    % Initialise Figure
    fig = fig + 1;
    
    switch fieldType

        case 'scalar'
            figName = [location, '_POD_', PODvar, '_Mode_Energy_Content'];

        case 'vector'
            figName = [location, '_POD_', cell2mat(PODvar), '_Mode_Energy_Content'];

    end
    
    set(figure(fig), 'outerPosition', [25, 25, 1275, 850], 'name', figName);
    set(gca, 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
             'fontSize', 20, 'layer', 'top');
    hold on;
    
    % Plot Modes
    cMap = viridis(3);
    bar(PODdata.modeEnergy, 0.75, ...
        'lineWidth', 2, 'faceColor', cMap(1,:));
    
    % Format Figure
    axis on;
    box on;
    grid off;
    xlim([0; 50]);
    ylim([0; 20]);
    tickData = (10:10:40);
    xticks(tickData);
    tickData = (4:4:16);
    yticks(tickData);
    xlabel('\bf{Mode}}', 'fontName', 'LM Roman 12');
    ylabel('{\bf{Energy Content (\it{%})}}', 'fontName', 'LM Roman 12');
    set(gca, 'outerPosition', [0.05, 0.05, 0.9, 0.9]);
    hold off;
    
    % Save Figure
    pause(2);
    
    exportgraphics(gcf, [userpath, '/Output/Figures/', figName, '.png'], 'resolution', 600);
    savefig(gcf, [userpath, '/Output/Figures/', figName, '.fig']);

end