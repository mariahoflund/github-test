OR_thrusterforce = 50e-3;   % N
SC_mass = 500;              % kg
thrust_acc = OR_thrusterforce/SC_mass;  % m/s^2
ext_ENV_muEarth = 3.986005e14; %earth gravitational constant, m^3/s^2

A = 42164500;           % m
psi_rads = 7.2921e-05;  % rad/s
psi_degday = 360.985647;    % deg/day

lambda_GEO = 79*pi/180;     % rad

%init orbit
sma0 = 42164200;        % m
% sma0 = 42170000;        % m
e0 = 0.0018;
inc0 = 0.1;             % deg
RAAN0 = 0;              % deg
Arg_of_per0 = 0;        % deg
nu0 = 203;              % deg

sim_start = juliandate(2020, 1, 1, 12, 0, 0);

[pos_start, vel_start] = keplerian2ijk(sma0,e0,inc0,RAAN0,Arg_of_per0,nu0);
