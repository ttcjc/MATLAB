%% PDA Reader v1.0

clear variables;
close all;
clc;

disp ('===============');
disp ('PDA Reader v1.0');
disp ('===============');
disp (' ');


%% Changelog

% v1.0 - Initial Commit


%% Do Stuff

load('~/MATLAB/CFD_Processing/particleData/PDA.mat');

massPerInjectorAnton = [1.996282029783190e-06;5.310740945893880e-06;1.128252209598640e-05;
                        2.838988653772200e-05;6.706500261665940e-05;1.515879037824680e-04;
                        3.438658384954620e-04;6.819114185828900e-04;0.001032379995671;0.001156210409242];

for i = 1:10
    injector = ['injector_', num2str(i)];
    injectors.(injector).distribution = injectorDistribution(PDA, 15 - i);
    injectors.(injector).thetaInner = injectorTheta(PDA, 16 - i);
    injectors.(injector).thetaOuter = injectorTheta(PDA, 15 - i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%

traverse = ['Traverse_', num2str(5)];
PDAdata = zeros(size(PDA.Pressure_110bar.Test2.(traverse).dataDroplet(:,7),1),4);
PDAdata(:,1) = PDA.Pressure_110bar.Test2.(traverse).dataDroplet(:,7); % Diameter [um]
PDAdata(:,2) = PDA.Pressure_110bar.Test2.(traverse).dataDroplet(:,5); % Axial Vel [m/s]
PDAdata(:,3) = PDA.Pressure_110bar.Test2.(traverse).dataDroplet(:,6); % Radial Vel [m/s]
PDAdata(:,4) = PDA.Pressure_110bar.Test2.(traverse).dataDroplet(:,4); % Residence Time [ms]

Umean = mean(PDAdata(:,2));
Vmean = mean(PDAdata(:,3));
velMeanA = sqrt(Umean ^ 2 + Vmean ^ 2)

velMag= sqrt(PDAdata(:,2) .^ 2 + PDAdata(:,3) .^ 2);
velMeanB = mean(velMag)
velWeightedMeanA = (sum(((4 / 3) * pi * ((PDAdata(:,1) * 1e-6 / 2) .^ 3)) .* velMag)) / sum(((4 / 3) * pi * ((PDAdata(:,1) * 1e-6 / 2) .^ 3)))

UweightedMean = (sum(((4 / 3) * pi * ((PDAdata(:,1) * 1e-6 / 2) .^ 3)) .* PDAdata(:,2))) / sum(((4 / 3) * pi * ((PDAdata(:,1) * 1e-6 / 2) .^ 3)));
VweightedMean = (sum(((4 / 3) * pi * ((PDAdata(:,1) * 1e-6 / 2) .^ 3)) .* PDAdata(:,3))) / sum(((4 / 3) * pi * ((PDAdata(:,1) * 1e-6 / 2) .^ 3)));
velWeightedMeanB = sqrt(UweightedMean ^ 2 + VweightedMean ^ 2)

PDAdata = sortrows(PDAdata,1);

i = 1;
while PDAdata(i,1) == 0
    PDAdata(i,:) = [];
end

PDAdata(:,1) = round(PDAdata(:,1));

i = 1;
while PDAdata(i,1) == 0
    PDAdata(i,1) = 1;
    i = i + 1;
end

velMag= sqrt(PDAdata(:,2) .^ 2 + PDAdata(:,3) .^ 2);
velWeightedMeanC = (sum(((4 / 3) * pi * ((PDAdata(:,1) * 1e-6 / 2) .^ 3)) .* velMag)) / sum(((4 / 3) * pi * ((PDAdata(:,1) * 1e-6 / 2) .^ 3)))

UweightedMean = (sum(((4 / 3) * pi * ((PDAdata(:,1) * 1e-6 / 2) .^ 3)) .* PDAdata(:,2))) / sum(((4 / 3) * pi * ((PDAdata(:,1) * 1e-6 / 2) .^ 3)));
VweightedMean = (sum(((4 / 3) * pi * ((PDAdata(:,1) * 1e-6 / 2) .^ 3)) .* PDAdata(:,3))) / sum(((4 / 3) * pi * ((PDAdata(:,1) * 1e-6 / 2) .^ 3)));
velWeightedMeanD = sqrt(UweightedMean ^ 2 + VweightedMean ^ 2)

classes = unique(PDAdata(:,1));
APL = zeros(size(classes,1),1);
VV = zeros(size(classes,1),1);

for i = 1:size(classes,1)
    index = find(PDAdata(:,1) == classes(i));
    APL(i) = (1 / sum(PDAdata(:,1) == classes(i))) * (sum((PDAdata(index,4) * 1e-3) .* velMag(index)));
end


%% Clean-up



%% Local Functions

function distribution = injectorDistribution(PDA, index)

    traverse = ['Traverse_', num2str(index)];
    PDAdiameters = sort(PDA.Pressure_110bar.Test4.(traverse).dataDroplet(:,7));
    
    i = 1;
    while PDAdiameters(i) == 0
        PDAdiameters(i) = [];
    end
    
    PDAdiameters = round(PDAdiameters);
    
    i = 1;
    while PDAdiameters(i) == 0
        PDAdiameters(i) = 1;
        i = i + 1;
    end
    
    distribution = zeros(size(unique(PDAdiameters),1),2);
    distribution(:,1) = unique(PDAdiameters);
    
    for i = 1:size(distribution,1)
        distribution(i,2) = (sum(PDAdiameters == distribution(i,1)) / size(PDAdiameters,1)) * 100;
    end
    
end

function theta = injectorTheta(PDA, index)

    traverse = ['Traverse_', num2str(index)];
    theta = atand(abs(PDA.Pressure_110bar.Test4.(traverse).Traverse / PDA.Pressure_110bar.Test4.(traverse).Zpos));
    
end