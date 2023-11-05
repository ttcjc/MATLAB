clA;


%%

radiusWheel = 300e-3;
widthWheel = 220e-3;
heightCentre = radiusWheel - 18e-3;

depthTread = 8e-3;
widthTreadA = 12e-3;
widthTreadB = 6e-3;
nTreadA = 2;
nTreadB = 2;

treadFillRatio = 1;

gamma = interp1([0, 10], [0.0757, 0.0742], 5);
rho = 1000;

U = 22.222222222222222;
w = U / radiusWheel;


%%

massTP = rho * U * ((nTreadA * (treadFillRatio * (widthTreadA * depthTread))) + (nTreadB * (treadFillRatio * (widthTreadB * depthTread))))

filmThickness = 0.1e-3; % filmThickness = (gamma / rho) * (1 / ((w * radiusWheel)^2));
filmArea = (widthWheel * filmThickness) + ...
           (nTreadA * (widthTreadA * depthTread)) - (nTreadA * ((widthTreadA - (2 * filmThickness)) * depthTread)) + ...
           (nTreadB * (widthTreadB * depthTread)) - (nTreadB * ((widthTreadB - (2 * filmThickness)) * depthTread));
massCA = rho * U * filmArea

ratio = massCA / (massTP + massCA)

disp(' ');
disp(' ');


%%

QS_parcelsPerSecond = 25e6 * 0.75
QS_Tstar = 1.044 / 40;
QS_parcelsPerTstar = QS_parcelsPerSecond / (1 / QS_Tstar)

FS_Tstar = (4 * 1.044) / 22.222222222222222;
FS_parcelsPerSecond = (QS_parcelsPerTstar * (1 / FS_Tstar)) / 4

disp(' ');
disp(' ');


%%

parcelsTP = ceil(FS_parcelsPerSecond * (1 - ratio))
parcelsCA = ceil(FS_parcelsPerSecond * (ratio))

disp(' ');
disp(' ');


%%

radiusInj = radiusWheel + 2e-3;

minTPangle = 270 + acosd(heightCentre / radiusInj) + 1
maxTPangle = minTPangle + 20

minCAangle = minTPangle
maxCAangle = (minCAangle - 1) + (360 - (2 * (acosd(heightCentre / radiusInj)))) - 1