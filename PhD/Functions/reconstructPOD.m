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

function [reconData] = reconstructPOD(reconData, PODvar, nModes, Nt, Ns, fieldType)
    
    disp('    Performing Field Reconstruction...');
    
    for i = nModes
        % Initialise Progress Bar
        wB = waitbar(0, ['Adding Mode #', num2str(i), ' to Reconstruction'], 'name', 'Progress');
        wB.Children.Title.Interpreter = 'none';
        
        % Identify Mode Contribution
        mode = ['M', num2str(i)];
        
        if strcmp(fieldType, 'scalar')
            reconData.(mode).modeMatrix = PODdata.A_coeff(:,i) * PODdata.phi_mode(:,i)';
            reconData.(mode).prime = cell(Nt,1);
            
            for j = 1:Nt
                reconData.(mode).prime{j} = zeros(Ns,1);
                
                for k = 1:Ns
                    reconData.(mode).prime{j}(k) = reconData.(mode).modeMatrix(j,k);
                end
                
                waitbar((j / Nt), wB);
            end
        
        elseif strcmp(fieldType, 'vector')
            uModeMatrix = PODdata.A_coeff(:,i) * PODdata.phi_mode((1:Ns),i)';
            vModeMatrix = PODdata.A_coeff(:,i) * PODdata.phi_mode(((Ns + 1):(2 * Ns)),i)';
            wModeMatrix = PODdata.A_coeff(:,i) * PODdata.phi_mode((((2 * Ns) + 1):end),i)';
            
            reconData.(mode).u.prime = cell(Nt,1);
            reconData.(mode).v.prime = reconData.(mode).u.prime;
            reconData.(mode).w.prime = reconData.(mode).u.prime;
            
            for j = 1:Nt
                reconData.(mode).u.prime{j} = zeros(Ns,1);
                reconData.(mode).v.prime{j} = reconData.(mode).u.prime{j};
                reconData.(mode).w.prime{j} = reconData.(mode).u.prime{j};
                
                for k = 1:Ns
                    reconData.(mode).u.prime{j}(k) = uModeMatrix(j,k);
                    reconData.(mode).v.prime{j}(k) = vModeMatrix(j,k);
                    reconData.(mode).w.prime{j}(k) = wModeMatrix(j,k);
                end
                
                 waitbar((j / Nt), wB);
            end

            reconData.(mode).modeMatrix = [uModeMatrix, vModeMatrix, wModeMatrix];
        else
            error('Unexpected Field Type');
        end
        
        delete(wB);
        
        % Add Mode to Reconstruction
        if strcmp(fieldType, 'scalar')
            
            for j = 1:Nt
                reconData.(PODvar).inst{j} = reconData.(PODvar).inst{j} + reconData.(mode).prime{j};
            end

            elseif strcmp(fieldType, 'vector')
                
                for j = 1:Nt
                    reconData.(PODvar{1}).inst{j} = reconData.(PODvar{1}).inst{j} + reconData.(mode).u.prime{j};
                    reconData.(PODvar{2}).inst{j} = reconData.(PODvar{2}).inst{j} + reconData.(mode).v.prime{j};
                    reconData.(PODvar{3}).inst{j} = reconData.(PODvar{3}).inst{j} + reconData.(mode).w.prime{j};
                end

        end
    
    end

end