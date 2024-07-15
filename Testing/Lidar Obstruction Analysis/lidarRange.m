%% Preamble

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

if exist('/mnt/Processing/Data', 'dir')
    saveLocation = '/mnt/Processing/Data';
else
    saveLocation = '~/Data';
end

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallelisation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures


%% Load Test Data

load('/mnt/Processing/Data/Numerical/MATLAB/volumeObstruction/Windsor_SB_fullScale_multiPhase/farWake/T1002_T3200_F50_D20_D400/Plane_YZ_19465_Origin_18637_0_076.mat');


%% Mie Efficiency Calculations

refIndex_ratio = 1.32352 + (5.15e-7)*i;
waveNumber = 905e-9 / 1.00028;

% % Visualise Broadband Efficiency Trends
% particleD = (logspace(-8, -2, 1000))';
% Q_e = zeros(height(particleD),1);
% Q_b = zeros(height(particleD),1);
% 
% for j = 1:height(particleD)
%     sizeParameter = (pi * particleD(j)) / waveNumber;
%     
%     result = mie(refIndex_ratio, sizeParameter);
%     
%     Q_e(j) = result(4);
%     Q_b(j) = result(7);
% end
% clear i;
% 
% % Initialise Figure
% fig = fig + 1;
% set(figure(fig), 'name', 'Mie Efficiencies', 'color', [1, 1, 1], ...
%                  'outerPosition', [25, 25, 650, 650], 'units', 'pixels')
% set(gca, 'positionConstraint', 'outerPosition', ...
%          'lineWidth', 2, 'fontName', 'LM Mono 12', 'fontSize', 16, 'layer', 'top', ...
%          'xScale', 'log', 'yScale', 'log');
% hold on;
% 
% % Plot Data
% plot(particleD, Q_e, 'b-', 'lineWidth', 1.5);
% plot(particleD, Q_b, 'r-', 'lineWidth', 1.5);
% xline(1e-6);
% xline(147e-6);
% 
% % Format Figure
% axis on;
% box on;
% grid off;
% xlim([1e-8; 1e-2]);
% ylim([1e-7; 1e2]);
% hold off;

% Calculate Relevant Efficiencies
D = (20e-6:2e-7:400e-6)';
Qe_range = zeros([height(D),1]);
Qb_range = zeros([height(D),1]);

for k = 1:height(D)
    sizeParameter = (pi * D(k)) / waveNumber;
    
    result = mie(refIndex_ratio, sizeParameter);
    
    Qe_range(k) = result(4);
    Qb_range(k) = result(7);
end
clear i;

QeInterp = griddedInterpolant(D, Qe_range);
QbInterp = griddedInterpolant(D, Qb_range);


%% Lidar Range Calculations

c = 299792458 / 1.00028;
alphaClearSky = 100;

R_max = 18;
% R_max = 50;

% dL = 0.1;

dT = (2 * (R_max / c)) / 100;

for i = 3 % [1,2,3] % 1:height(obstructData.raysOfInterest)
    
    maxPing = zeros([51,1]);
    
    for j = 1:10 % 1:height(obstructData.time)
%         R_o = height(obstructData.inst.rayData.samplePoints{j}{i}) * dL;
%         R_o = 40;
        R_o = 16;
        A_o = inf;

%         axisDisp = 22.25e-3;
        axisDisp = 27e-3;
        D_t0 = 4e-3;
        D_r0 = 40e-3;
        gamma_t = deg2rad(0.06);
        gamma_r = deg2rad(0.1);
        A_r = (tau * D_r0^2) / 4;
        eta_t = 1;
        eta_r = 1;
        Ppeak = 1;
%         pulseWidth = 2.5e-10;
        pulseWidth = 5e-9;

        R1 = (axisDisp - (D_t0 / 2) - (D_r0 / 2)) / (tan(gamma_t / 2) + tan(gamma_r / 2));
        R2 = (axisDisp - (D_r0 / 2) + (D_t0 / 2)) / (tan(gamma_r / 2) - tan(gamma_t / 2));

        Gamma_o = 0.2;

        sampleDist = (dL:dL:R_max)';
        alpha = zeros([height(sampleDist),1]);
        beta = alpha;
        H_t = alpha;
        H_c = alpha;
        H = alpha;
        Phi_t = alpha;
        Phi_r = alpha;
        zeta = alpha;
        P_r = alpha;

        nParticles = zeros([height(sampleDist),1]); nParticles(1:height(obstructData.inst.rayData.nParticles{j}{i})) = obstructData.inst.rayData.nParticles{j}{i}; % * 60;
        dParticles = zeros([height(sampleDist),1]); dParticles(1:height(obstructData.inst.rayData.d{j}{i})) = obstructData.inst.rayData.d{j}{i} * 1e-6; % * 5e-6;
        massParticles = zeros([height(sampleDist),1]); massParticles(1:height(obstructData.inst.rayData.density{j}{i})) = obstructData.inst.rayData.density{j}{i};

        for k = 1:height(sampleDist)
            R = sampleDist(k);

            D_t = D_t0 + (2 * R * tan(gamma_t / 2));
            D_r = D_r0 + (2 * R * tan(gamma_r / 2));
            beamA = (tau * D_t^2) / 4;

            D1 = D_t0 + (2 * (R - (dL / 2)) * tan(gamma_t / 2));
            D2 = D_t0 + (2 * (R + (dL / 2)) * tan(gamma_t / 2));
            beamV = (tau / 6) * dL * ((D1^2 / 4) + ((D1 * D2) / 4) + (D2^2 / 4));

            nParticles(k) = nParticles(k) * (beamV / cellSize.volume);
            massParticles(k) = massParticles(k) * beamV;
            
            Qe = QeInterp(dParticles(k));
            Qb = QbInterp(dParticles(k));
            
            sigma_e = Qe * ((tau * (dParticles(k))^2) / 8);
            sigma_b = Qb * ((tau * (dParticles(k))^2) / 8);
            
            alpha(k) = alphaClearSky + ((nParticles(k) * sigma_e) / beamV);

            beta(k) = (nParticles(k) * sigma_b) / beamV;

            if A_o >= beamA
                H_t(k) = (Gamma_o * (diracDelta((R - R_o), 8) / dT)) + beta(k);
            else
                H_t(k) = (Gamma_o * (diracDelta((R - R_o), 8) / dT) * (At / beamA)) + beta(k);
            end

            Phi_t(k) = 2 * acos(((D_t / 2)^2 - (D_r / 2)^2 + axisDisp^2) / (2 * axisDisp * (D_t / 2)));
            Phi_r(k) = 2 * acos(((D_r / 2)^2 - (D_t / 2)^2 + axisDisp^2) / (2 * axisDisp * (D_r / 2)));

            if R <= R1
                zeta(k) = 0;
            elseif R >= R2
                zeta(k) = 1;
            else
                zeta(k) = (((D_t / 2)^2 * (Phi_t(k) - sin(Phi_t(k)))) + ((D_r / 2)^2 * (Phi_r(k) - sin(Phi_r(k))))) / (tau * (D_t / 2)^2);
            end

            H_c(k) = exp(-sum(alpha(1:k) * dL))^2 * (zeta(k) / (tau * R^2));

            H(k) = H_c(k) + H_t(k);

            if k == 1
                Hinterp = griddedInterpolant([0; sampleDist(1)], [0; H(1)]);
            else
                Hinterp = griddedInterpolant(sampleDist(1:k), H(1:k));
            end

            t = (0:dT:(2 * (R / c)))';
            
            Pt_conv = zeros([height(t),1]);
            H_conv = Pt_conv;
            
            for ii = 1:height(t)

                if t(ii) <= pulseWidth
                    Pt_conv(ii) = Ppeak * (sin((tau / 2) * (t(ii) / pulseWidth)))^2;
                end

                H_conv(ii) = Hinterp(R - (c * (t(ii) / 2)));
            end
            clear ii;

            convolution = sum(Pt_conv .* H_conv * dT);
            
            P_r(k) = (A_r * eta_t * eta_r) * convolution;
        end
        clear k;
        
        index= find(sampleDist == R_o);
        
%         P_r = P_r / max(P_r); P_r = P_r(1:index); P_r(P_r < 1e-14) = 0; P_r((end-2):end) = 0;
        sampleDist = sampleDist / R_o; % sampleDist = sampleDist(1:index);
        nParticles = nParticles(1:index);
        dParticles = dParticles(1:index);
        
        figTime = num2str(obstructData.time(j), ['%.', num2str(timePrecision), 'f']);
        
        figTitle = '-';
        figSubtitle = [figTime, ' \it{s}'];
        
        %%%%
        
        % Initialise Figure
        fig = fig + 1;
        figName = ['Lidar_Detections_Ray_', obstructData.raysOfInterest{i}, '_T', erase(figTime, '.')];
        set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                         'outerPosition', [25, 25, 650, 650], 'units', 'pixels')
        set(gca, 'positionConstraint', 'outerPosition', ...
                 'lineWidth', 2, 'fontName', 'LM Mono 12', 'fontSize', 16, 'layer', 'top');
        hold on;

        % Plot Data
        plot(sampleDist, P_r, 'lineWidth', 1.5, 'color', ([74, 24, 99] / 255));
        xline(1, '--', 'Intended Target', 'labelVerticalAlignment', 'middle', ...
                                               'lineWidth', 1.5, 'color', [0, 0, 0]);
        
        % Format Figure
%         title(figTitle);
%         subtitle(figSubtitle);
        box on;
        grid off;
        hold off;
%         xlim([-0.1; 1.1]);
%         ylim([-7.2e-9; 7.92e-8]);
%         xlim([-0.1; 1.1]);
%         ylim([-0.1; 1.1]);
%         xticks(0:0.25:1);
%         yticks(0:1.2e-8:7.2e-8);
%         xticks(0:0.25:1);
%         yticks(0:0.25:1);
        xticks([]);
        yticks([]);
%         xtickformat('%0.2f');
%         ytickformat('%0.1f');
%         xtickformat('%0.2f');
%         ytickformat('%0.2f');
%         xlabel({'{\bf{Normalised Distance}}'; '-'}, 'fontName', 'LM Roman 12');
%         ylabel({'-'; '{\bf{Relative Return Signal}}'}, 'fontName', 'LM Roman 12');
        xlabel({'{\bf{Distance}}'}, 'fontName', 'LM Roman 12');
        ylabel({'{\bf{Return Signal}}'}, 'fontName', 'LM Roman 12');

        % Save Figure
        pause(1);
        exportgraphics(gcf, [userpath, '/Output/Figures/', figName, '.png'], 'resolution', 600);
        close(fig);
        
        %%%%
        
%         % Initialise Figure
%         fig = fig + 1;
%         figName = ['Particle_Count_Ray_', obstructData.raysOfInterest{i}, '_T', erase(figTime, '.')];
%         set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
%                          'outerPosition', [25, 25, 650, 650], 'units', 'pixels')
%         set(gca, 'positionConstraint', 'outerPosition', ...
%                  'lineWidth', 2, 'fontName', 'LM Mono 12', 'fontSize', 16, 'layer', 'top');
%         hold on;
% 
%         % Plot Data
%         plot(sampleDist, nParticles, 'b-', 'lineWidth', 1.5, 'color', ([74, 24, 99] / 255));
%         xline(1, '--', 'Intended Target', 'labelVerticalAlignment', 'middle', ...
%                                                'lineWidth', 1.5, 'color', [0, 0, 0]);
%         
%         % Format Figure
%         title(figTitle);
%         subtitle(figSubtitle);
%         box on;
%         grid off;
%         hold off;
%         xlim([-0.1; 1.1]);
%         ylim([-40; 440]);
%         xticks(0:0.25:1);
%         yticks(0:80:400);
%         xtickformat('%0.2f');
%         ytickformat('%0.f');
%         xlabel({'{\bf{Normalised Distance}}'; '-'}, 'fontName', 'LM Roman 12');
%         ylabel({'-'; '{\bf{Particle Count}}'}, 'fontName', 'LM Roman 12');
% 
%         % Save Figure
%         pause(1);
%         exportgraphics(gcf, [userpath, '/Output/Figures/', figName, '.png'], 'resolution', 600);
%         close(fig);

        %%%%
        
%         % Initialise Figure
%         fig = fig + 1;
%         figName = ['Particle_Diameters_Ray_', obstructData.raysOfInterest{i}, '_T', erase(figTime, '.')];
%         set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
%                          'outerPosition', [25, 25, 650, 650], 'units', 'pixels')
%         set(gca, 'positionConstraint', 'outerPosition', ...
%                  'lineWidth', 2, 'fontName', 'LM Mono 12', 'fontSize', 16, 'layer', 'top');
%         hold on;
% 
%         % Plot Data
%         plot(sampleDist, (dParticles * 1e6), 'b-', 'lineWidth', 1.5, 'color', ([74, 24, 99] / 255));
%         xline(1, '--', 'Intended Target', 'labelVerticalAlignment', 'middle', ...
%                                                'lineWidth', 1.5, 'color', [0, 0, 0]);
%         
%         % Format Figure
%         title(figTitle);
%         subtitle(figSubtitle);
%         box on;
%         grid off;
%         hold off;
%         xlim([-0.1; 1.1]);
%         ylim([-15; 165]);
%         xticks(0:0.25:1);
%         yticks(0:25:150);
%         xtickformat('%0.2f');
%         ytickformat('%0.f');
%         xlabel({'{\bf{Normalised Distance}}'; '-'}, 'fontName', 'LM Roman 12');
%         ylabel({'-'; '{\bf{Mean Particle Diameter (\mum)}}'}, 'fontName', 'LM Roman 12');
% 
%         % Save Figure
%         pause(1);
%         exportgraphics(gcf, [userpath, '/Output/Figures/', figName, '.png'], 'resolution', 600);
%         close(fig);
    end
    clear j;
    
end
clear i;