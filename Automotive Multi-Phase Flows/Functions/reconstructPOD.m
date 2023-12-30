%% Snapshot POD Field Reconstruction Tool v2.0
% ----
% Reconstructs Field Using POD Data Produced Using 'performPOD'
%
% J. Weiss
% "A Tutorial on the Proper Orthogonal Decomposition"
% 2019 AIAA Aviation Forum, 17-21 June 2019, Dallas, Texas, United States
% ----
% Usage: reconData = reconstructPOD(reconData, PODdata, field, nModes, Ns, Nt, fieldType, saveModes)
%
%        'reconData' -> Structure Containing Position and Field Data
%        'PODdata'   -> Structure Containing Previously Processed POD Data
%        'field'     -> Field Variable Used to Perform POD Stored as String
%        'Ns'        -> Number of Spatial Points
%        'Nt'        -> Number of Temporal Points
%        'fieldType' -> Desired Field Type Stored as String
%        'saveModes' -> Save Individual Mode Contributions [True/False]


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Formatting Updates
% v1.2 - Further Formatting Updates
% v2.0 - Update To Support Format of Separated Reconstruction Scripts


%% Supported Field Types

% Fluctuating Scalar Field: 'scalar'
% Fluctuating Vector Field: 'vector'


%% Main Function

function reconData = reconstructPOD(reconData, field, nModes, Ns, Nt, fieldType, saveModes)
    
    disp(['    Performing Field Reconstruction Using ' num2str(width(nModes)), ' Mode(s)...']);
    
    % Initialise Progress Bar
    wB = waitbar(0, 'Adding Mode(s) to Reconstruction', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    
    for i = nModes        
        mode = ['M', num2str(i)];
        
        % Identify Mode Contribution        
        switch fieldType

            case 'scalar'
                modeMatrix = reconData.POD.alpha(:,i) * reconData.POD.phi(:,i)';
                reconData.POD.(mode).prime = cell(Nt,1);
                
                for j = 1:Nt
                    reconData.POD.(mode).prime{j} = modeMatrix(j,:)';
                end
                clear j;

            case 'vector'
                uModeMatrix = reconData.POD.alpha(:,i) * reconData.POD.phi((1:Ns),i)';
                vModeMatrix = reconData.POD.alpha(:,i) * reconData.POD.phi(((Ns + 1):(2 * Ns)),i)';
                wModeMatrix = reconData.POD.alpha(:,i) * reconData.POD.phi((((2 * Ns) + 1):end),i)';
                
                reconData.POD.(mode).u.prime = cell(Nt,1);
                reconData.POD.(mode).v.prime = reconData.(mode).u.prime;
                reconData.POD.(mode).w.prime = reconData.(mode).u.prime;
                
                for j = 1:Nt
                    reconData.POD.(mode).u.prime{j} = uModeMatrix(j,:)';
                    reconData.POD.(mode).v.prime{j} = vModeMatrix(j,:)';
                    reconData.POD.(mode).w.prime{j} = wModeMatrix(j,:)';
                end
                clear j;
        
        end
        
        % Add Mode to Reconstruction
        switch fieldType

            case 'scalar'
                
                for j = 1:Nt
                    reconData.(field).recon.inst{j} = reconData.(field).recon.inst{j} + ...
                                                                reconData.POD.(mode).prime{j};
                end

            case 'vector'
                
                for j = 1:Nt
                    reconData.(field{1}).recon.inst{j} = reconData.(field{1}).recon.inst{j} + ...
                                                          reconData.POD.(mode).u.prime{j};
                    reconData.(field{2}).recon.inst{j} = reconData.(field{2}).recon.inst{j} + ...
                                                          reconData.POD.(mode).v.prime{j};
                    reconData.(field{3}).recon.inst{j} = reconData.(field{3}).recon.inst{j} + ...
                                                          reconData.POD.(mode).w.prime{j};
                end
                clear j;
        
        end
        
        if ~saveModes
            reconData.POD = rmfield(reconData.POD, mode);
        end
        
        % Update Waitbar
        waitbar((i / length(nModes)), wB);
    end
    clear i;
    
    delete(wB);

end