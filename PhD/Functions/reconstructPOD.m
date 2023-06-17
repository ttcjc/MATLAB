%% Snapshot POD Field Reconstruction Tool v1.1
% ----
% Reconstructs Field Using POD Data Produced Using 'performPOD'
%
% J. Weiss
% "A Tutorial on the Proper Orthogonal Decomposition"
% 2019 AIAA Aviation Forum, 17-21 June 2019, Dallas, Texas, United States
% ----
% Usage: reconData = reconstructPOD(reconData, PODdata, PODvar, nModes, Ns, Nt, fieldType, saveModes)
%
%        'reconData' -> Structure Containing Position and Field Data
%        'PODdata'   -> Structure Containing Previously Processed POD Data
%        'PODvar'    -> Field Variable Used to Perform POD Stored as String
%        'Ns'        -> Number of Spatial Points
%        'Nt'        -> Number of Temporal Points
%        'fieldType' -> Desired Field Type Stored as String
%        'saveModes' -> Save Individual Mode Contributions [True/False]


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Formatting Updates


%% Supported Field Types

% Fluctuating Scalar Field: 'scalar'
% Fluctuating Vector Field: 'vector'


%% Main Function

function reconData = reconstructPOD(reconData, PODdata, PODvar, nModes, Ns, Nt, fieldType, saveModes)
    
    disp(['    Performing Field Reconstruction Using ' num2str(width(nModes)), ' Mode(s)...']);
    
    % Initialise Progress Bar
    wB = waitbar(0, 'Adding Mode(s) to Reconstruction', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    
    for i = nModes        
        mode = ['M', num2str(i)];
        
        % Identify Mode Contribution        
        switch fieldType

            case 'scalar'
                modeMatrix = PODdata.A_coeff(:,i) * PODdata.phi_mode(:,i)';
                reconData.(mode).prime = cell(Nt,1);
                
                for j = 1:Nt
                    reconData.(mode).prime{j} = modeMatrix(j,:)';
                end
                clear j;

            case 'vector'
                uModeMatrix = PODdata.A_coeff(:,i) * PODdata.phi_mode((1:Ns),i)';
                vModeMatrix = PODdata.A_coeff(:,i) * PODdata.phi_mode(((Ns + 1):(2 * Ns)),i)';
                wModeMatrix = PODdata.A_coeff(:,i) * PODdata.phi_mode((((2 * Ns) + 1):end),i)';
                
                reconData.(mode).u.prime = cell(Nt,1);
                reconData.(mode).v.prime = reconData.(mode).u.prime;
                reconData.(mode).w.prime = reconData.(mode).u.prime;
                
                for j = 1:Nt
                    reconData.(mode).u.prime{j} = uModeMatrix(j,:)';
                    reconData.(mode).v.prime{j} = vModeMatrix(j,:)';
                    reconData.(mode).w.prime{j} = wModeMatrix(j,:)';
                end
                clear j;
        
        end
        
        % Add Mode to Reconstruction
        switch fieldType

            case 'scalar'
                
                for j = 1:Nt
                    reconData.(PODvar).inst{j} = reconData.(PODvar).inst{j} + reconData.(mode).prime{j};
                end

            case 'vector'
                
                for j = 1:Nt
                    reconData.(PODvar{1}).inst{j} = reconData.(PODvar{1}).inst{j} + reconData.(mode).u.prime{j};
                    reconData.(PODvar{2}).inst{j} = reconData.(PODvar{2}).inst{j} + reconData.(mode).v.prime{j};
                    reconData.(PODvar{3}).inst{j} = reconData.(PODvar{3}).inst{j} + reconData.(mode).w.prime{j};
                end
                clear j;
        
        end
        
        if ~saveModes
            reconData = rmfield(reconData, mode);
        end
        
        % Update Waitbar
        waitbar((i / length(nModes)), wB);
    end
    clear i;
    
    delete(wB);

end