function [M] = DyneinAndNexinContribution(obj, i)

% Initializing parameters using values in [35], [36]

a = 200e-9; % Axoneme diameter.

rho = 10e-3; % Mean number density of motors.

pi0p = 0.17; % Characteristic timescale for attachment (for n+).

pi0m = 0.25; % Characteristic timescale for attachment (for n-).

eps0p = 0.73; % Characteristic timescale for detachment (for n+).

eps0m = 0.75; % Characteristic timescale for detachment (for n-).

fc = 1e-12; % Charactertistic (unbinding) force, about which motors detach exponentially fast.
    % Should be in range (0.5,2.5)pN.

f0 = 2e-12; % Stall force of dynein motor.
    % Should be in range (1-5)pN.

v0 = 5e-6;  % Motor walking speed at zero load.
    % Should be in range (5-7) micrometers.

K = 2e3;  % Stiffness of Nexin links (in Chlamydomonas). 

phi0 = pi/2; % Value of the angle, phi, the centerline makes with the x-axis at s=0


% Initializing variables 
Qtemp = obj.Q;
N = obj.Np;
dL = obj.DeltaL;
L = N * dL; % Length of filament.

%% Step 2
t = QuaternionRotation(Qtemp(i,1:4),[1;0;0]);
phi = acos(t(1));
sliding_displacement = a*(phi - phi0);

%% Step 3
dphids = 
dsliding_displacementdt = 

end