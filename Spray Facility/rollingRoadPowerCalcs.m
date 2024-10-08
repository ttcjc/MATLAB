run preamble;


%%

g = 9.81; % m / s^2


%%

mass_model = 325; % kg (25% 2023 Ford Fiesta Kerb Weight)
mu_platen = 1; % (Worst Case Scenario)
velocity_road = 35; % m / s
radius_driven = 0.0762; % m
velocity_motor = 29.16; % rev / s
efficiency_drive = 0.7;
inertia_drive = 1.28; % kg / m^2
time_ramp = 60; % s
safetyFactor = 1.1;


%% 

F_gravity = mass_model * g;
F_suction = 0;
F_lift = 0;

F_vertical = F_gravity + F_suction - F_lift;


%%

F_friction = mu_platen * F_vertical;

F_horizontal = F_friction;


%%

velocity_Ratio = velocity_road / (tau * radius_driven * velocity_motor);


%%

torque = (((abs(F_horizontal) * (velocity_Ratio * radius_driven)) / efficiency_drive) + ...
         (inertia_drive * ((velocity_road / radius_driven) / time_ramp))) * safetyFactor
     
power = torque * ((velocity_road / radius_driven) / velocity_Ratio)