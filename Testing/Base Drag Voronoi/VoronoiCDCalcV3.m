function [AverageCPVoronoi,FigNo] = VoronoiCDCalcV3(TappingCoords,BoundingBox,PressureData,FigNo,Range,FrontOrBase,PlotReq,plotType,tileSelect)
% VoronoiCDCalc Calculates Pressure CD using the voronoi method
%   Inputs:
%   - TappingCoords - Column 1 is X, Column 2 is Y
%   - BoundingBox - Outer limits for voronoi plot
%   - PressureData - CP Values should be in column 1
%   - FigNo - Integer for figure number
%   - Range - [LowerLim UpperLim] or 'auto'
%   - FrontOrBase - 0 is front, 1 is base, used for inverting x coordinates
%   on the base
%   - plot required 1 = on, 0 = off
%   - plot type ('main') for one single plot, ('sub') is multiple voronois
%   on one plot needed
% V3 change adds subplot

warning off
[V,C,XY] = VoronoiLimitJamesBodge(TappingCoords(:,1), TappingCoords(:,2), 'bs_ext', BoundingBox, 'figure', 'off');
warning on
[LIA, LOCB] = ismember(TappingCoords,XY,'rows');

LOCB = LOCB(LOCB ~= 0);
PressureTapsSorted = XY(LOCB,:);



C_Sorted = C(LOCB);

if strcmp(FrontOrBase,'front')
    V(:,1) = V(:,1);
else
    V(:,1) = -V(:,1);
end
dA = zeros(height(TappingCoords),1);
for i = 1:height(TappingCoords)
    dA(i) = polyarea(V(C_Sorted{i},1),V(C_Sorted{i},2));
end
TotalArea = sum(dA);


AreaWeightedPressure = zeros(length(PressureData),1);
for i = 1:length(PressureData)
    AreaWeightedPressure(i) = PressureData(i,1)*dA(i);
end
AverageCPVoronoi = sum(AreaWeightedPressure)/TotalArea;

%% Plotting pressure in each region of voronoi diagram


if PlotReq == 1
    if strcmp(plotType,'main')
        figure(FigNo)
        FigNo = FigNo+1;
        hold on
        for i = 1:height(C_Sorted)
            x = V(C_Sorted{i},1);
            y = V(C_Sorted{i},2);
            a = PressureData(i,1);
            patch(x,y,a);
        
        end
        colormap('parula');
        caxis(Range);
        
        c = colorbar;
        c.Label.String = 'Cp';
        xlabel('y (m)');
        ylabel('z (m)');
        ytickformat('%.2f');
        xtickformat('%.2f');
        %         hold on
        %         plot(boundingBox(:,1), boundingBox(:,2), 'color', 'r', 'lineWidth', 2);
        scatter(TappingCoords(:,1), TappingCoords(:,2), 25, 'k', 'filled');
        set(gca, 'DataAspectRatio', [1,1,1]);
    elseif strcmp(plotType,'sub')
        figure(FigNo)
        FigNo = FigNo+1;
        if tileSelect == 1
            subplot(1,2,1)
        else
            subplot(1,2,2)
        end
        hold on
        for i = 1:height(C_Sorted)
            x = V(C_Sorted{i},1);
            y = V(C_Sorted{i},2);
            a = PressureData(i,1);
            patch(x,y,a);
        
        end
        colormap('parula');
        caxis(Range);
        
        c = colorbar;
        c.Label.String = 'Cp';
        xlabel('y (m)');
        ylabel('z (m)');
        ytickformat('%.2f');
        xtickformat('%.2f');
        %         hold on
        %         plot(boundingBox(:,1), boundingBox(:,2), 'color', 'r', 'lineWidth', 2);
        scatter(TappingCoords(:,1), TappingCoords(:,2), 25, 'k', 'filled');
        set(gca, 'DataAspectRatio', [1,1,1]);
    else
    end

else
end



end
