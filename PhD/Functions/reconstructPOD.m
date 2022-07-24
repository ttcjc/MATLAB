%% Snapshot POD Field Reconstruction Tool v1.0
% ----
% Reconstructs Field Using POD Data Produced Using 'performPOD'
%
% Weiss, Julien: A Tutorial on the Proper Orthogonal Decomposition. In: 2019 AIAA Aviation Forum. 17–21
% June 2019, Dallas, Texas, United States.
% ----
% Usage: [] = reconstructPOD()


%% Changelog

% v1.0 - Initial Commit


%% Supported Field Types

% Fluctuating Scalar Field: 'scalar'
% Fluctuating Vector Field: 'vector'


%% Main Function

function reconData = reconstructPOD(reconData, PODdata, PODvar, nModes, Ns, Nt, fieldType)
    
    disp('    Identifying Modal Contributions to Fluctuating Field...');
    
    for i = nModes
        % Initialise Progress Bar
        wB = waitbar(0, ['Identifying Mode #', num2str(i), ' Contributions'], 'name', 'Progress');
        wB.Children.Title.Interpreter = 'none';
        
        % Identify Mode Contribution
        mode = ['M', num2str(i)];
        
        switch fieldType

            case 'scalar'
                modeMatrix = PODdata.A_coeff(:,i) * PODdata.phi_mode(:,i)';
                reconData.(mode).prime = cell(Nt,1);
                
                for j = 1:Nt
                    reconData.(mode).prime{j} = modeMatrix(j,:)';
                    
                    waitbar((j / Nt), wB);
                end

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
                    
                    waitbar((j / Nt), wB);
                end
                
                modeMatrix = [uModeMatrix, vModeMatrix, wModeMatrix]; %#ok<NASGU>
            
            otherwise
                error('Unrecognised Field Type');
        
        end
        
        delete(wB);
    end

    disp(' ');

    disp('    Performing Field Reconstruction...');

    for i = nModes
        % Initialise Progress Bar
        wB = waitbar(0, ['Adding Mode #', num2str(i), ' to Reconstruction'], 'name', 'Progress');
        wB.Children.Title.Interpreter = 'none';
        
        % Add Mode to Reconstruction
        switch fieldType

            case 'scalar'
                
                for j = 1:Nt
                    reconData.(PODvar).inst{j} = reconData.(PODvar).inst{j} + reconData.(mode).prime{j};

                    waitbar((j / Nt), wB);
                end

            case 'vector'
                
                for j = 1:Nt
                    reconData.(PODvar{1}).inst{j} = reconData.(PODvar{1}).inst{j} + reconData.(mode).u.prime{j};
                    reconData.(PODvar{2}).inst{j} = reconData.(PODvar{2}).inst{j} + reconData.(mode).v.prime{j};
                    reconData.(PODvar{3}).inst{j} = reconData.(PODvar{3}).inst{j} + reconData.(mode).w.prime{j};

                    waitbar((j / Nt), wB);
                end
        
        end

        delete(wB);
    end

end