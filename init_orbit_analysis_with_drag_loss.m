clear all
%close all
% sim_SC_init_MoI_EOR_wet=[
%     1789.4664 -7.6657 14.0713
%     -7.6657 815.5687 -2.1297
%     14.0713 -2.1297 1862.6698];

% sim_SC_init_MoI_EOR_wet=[
%     1736.3548 -8.8513 -4.0280
%     -8.8513 567.9422 -2.0455
%     -4.0280 -2.0455 1668.1628];

sim_SC_init_MoI_EOR_wet=[
    1600 0 0;
    0 500 0;
    0 0 1500];

% COM = [0, 0, 468.85e-3];
% COM = [0.063, -0.004, 0.25+0.3];   %Inmarsat COM plus margin
COM = [0, 0, 0.25+0.3];



perigee_vector_eci = [0, 1, 0]; %update later


ext_ENV_muEarth = 3.986005e14; %earth gravitational constant, m^3/s^2
P_solar = 4.56e-6;  % solar radiation pressure
R_E = 6.3781e6;     %Earth radius

% Attitude
attitude_string = "z nadir, SA towards sun";
sun_pointing = 1;   % 0: sun pointing out of plane, 1: sun pointing in orbit plane
sun_pointing_SA = 0;    %0: sun ponting SA, 1: pointing in z+
theta = deg2rad(0);    %add rotation around z to attitude

% 3-body perturbation
enable_3body = "off";

% Atmospheric drag
drag_state = "off";     %with or without drag?
atmo_case = "average";      %average or max
enable = 1; %enable atmospheric model
alpha_acc = sqrt(0.7);
%geomagnetic activity
% aph_33_12 = [12, 6, 5, 9, 9, 6, 5, 9];
% aph_57_36 = [9, 5, 3, 12, 5, 6, 18];
% aph = [12, 18, 22, 22, 22, mean(aph_33_12), mean(aph_57_36)];


mission.StartDate = datetime(2020, 12, 21);     %2015-03-20 new mooon and equinox %2016-06-20 full moon
mission.StartDateJD = juliandate(mission.StartDate);
start_year = year(mission.StartDate);
start_doy = floor(mission.StartDateJD - juliandate(start_year, 1, 1))+1;


orbit = "SSTO";


switch(orbit)

    case "GTO"
        r_a = 3.5249e7+R_E;      %apogeumhöjd
        r_p = 2.5e5+R_E;    %perigeumhöjd
        a = (r_a + r_p)/2;  %semi major axis
        e = (r_a - r_p)/(r_a + r_p);      %eccentricity
        mission.Satellite.Inclination    = 6;    % deg
        mission.Satellite.ArgOfPeriapsis = 180;  % (deg) determines the orientation of inclination with respect to perigee
        mission.Satellite.RAAN           = 270;    % (deg) angle from vernal equinox to the ascending node (rotates the orbit around the z-axis)
        mission.Satellite.TrueAnomaly    = 0;    % deg

    case "SSTO"
        r_a = 6e7+R_E;      %apogeumhöjd
        r_p = 1.85e5+R_E;    %perigeumhöjd
        a = (r_a + r_p)/2;  %semi major axis
        e = (r_a - r_p)/(r_a + r_p);      %eccentricity
        mission.Satellite.Inclination    = 28.5;    % deg
        mission.Satellite.ArgOfPeriapsis = 180;  % (deg) determines the orientation of inclination with respect to perigee
        mission.Satellite.RAAN           = 270;    % (deg) angle from vernal equinox to the ascending node (rotates the orbit around the z-axis)
        mission.Satellite.TrueAnomaly    = 0;    % deg

    case "GEO"
        a = 42164e3;
        e = 0;

        % GEO longitude
        geo_lla = [0, 64, 0];  %lat, long, alt
        start_position = lla2eci(geo_lla, [year(mission.StartDate), month(mission.StartDate), day(mission.StartDate), ...
            hour(mission.StartDate), minute(mission.StartDate), second(mission.StartDate)]);
        start_angle = rad2deg(angle(start_position(2)*1i+start_position(1)));

        mission.Satellite.Inclination    = 0;    % deg
        mission.Satellite.ArgOfPeriapsis = 0;
        mission.Satellite.RAAN           = 0;
        mission.Satellite.TrueAnomaly    = start_angle;    % deg

end

switch(atmo_case)

    case "average"
        aph = 17*ones(1, 7);
        F107 = 150;

    case "max"
        aph = 200*ones(1,7);
        F107 = 300;

end
       

mission.Satellite.mass = 1025.94474+15;    %wet mass 
Cd = 2.179;     %drag coefficient
G = 6.6743e-11; %m3 kg-1 s-2

%Keplerian orbit period
T = 2*pi*sqrt(a^3/ext_ENV_muEarth);

mission.Satellite.area = 16;
mission.Duration = round(3*T);
% mission.Duration = 60*60*24*3;
% mission.Duration = 60*60*120;             
sim_step_size = 10;

mission.Satellite.SemiMajorAxis  = a;
mission.Satellite.Eccentricity   = e;

%% OPEN SIMULINK

mission.mdl = "orbit_propagator_matlab2023b";
open_system(mission.mdl);
mission.Satellite.blk = mission.mdl + "/Orbit Propagator";

%% SET PARAMETERS

set_param(mission.Satellite.blk, ...
    "startDate",      num2str(juliandate(mission.StartDate)), ...
    "stateFormatNum", "Orbital elements", ...
    "propagator",     "Numerical (high precision)", ...
    "orbitType",      "Keplerian", ...
    "semiMajorAxis",  "mission.Satellite.SemiMajorAxis", ...
    "eccentricity",   "mission.Satellite.Eccentricity", ...
    "inclination",    "mission.Satellite.Inclination", ...
    "raan",           "mission.Satellite.RAAN", ...
    "argPeriapsis",   "mission.Satellite.ArgOfPeriapsis", ...
    "trueAnomaly",    "mission.Satellite.TrueAnomaly");

set_param(mission.Satellite.blk, ...
    "centralBody",  "Earth", ...
    "outportFrame", "ICRF", ...
    "mass", "mission.Satellite.mass", ...
    "useDrag", drag_state, ...
    "useThirdBodyGravity", enable_3body ...
    );

set_param(mission.mdl, ...
    "SolverType", "Fixed-step", ...
    "FixedStep", "sim_step_size", ...
    "StopTime",   string(mission.Duration));

%%
mission.SimOutput = sim(mission.mdl);
%%
% mission.Satellite.TimeseriesPosICRF = mission.SimOutput.yout{1}.Values;
% mission.Satellite.TimeseriesVelICRF = mission.SimOutput.yout{2}.Values;
% 
% mission.Satellite.TimeseriesPosICRF.TimeInfo.StartDate = mission.StartDate;
% mission.Satellite.TimeseriesVelICRF.TimeInfo.StartDate = mission.StartDate;
% coe = mission.SimOutput.yout{10}.Values.Data;
% idx = 1;
% scenario = satelliteScenario;
% sat = satellite(scenario, mission.Satellite.TimeseriesPosICRF, mission.Satellite.TimeseriesVelICRF, ...
%     "CoordinateFrame", "inertial");
% 
%     % Retrieve states in geographic coordinates
%     [llaData, ~, llaTimeStamps] = states(sat(idx), "CoordinateFrame","geographic");
%     % Organize state data for each satellite in a seperate timetable
%     mission.Satellite.LLATable{idx} = timetable(llaTimeStamps', llaData(1,:)', llaData(2,:)', llaData(3,:)',...
%         'VariableNames', {'Lat_deg','Lon_deg', 'Alt_m'});
%     tout = mission.SimOutput.tout;
%     yout = mission.SimOutput.yout;



%%

var_string = [];

for i = 1:numElements(mission.SimOutput.yout)

eval(['var',num2str(i),' = mission.SimOutput.yout{i}.Values.Data;']);

time = mission.SimOutput.tout;
var_string = [var_string, string(mission.SimOutput.yout{i}.Name)];

end

% filename = "1rev " + orbit + " att " + attitude_string + " sun azim 105 alt " + num2str(sun_alt_vec(alt));
% save("angular momentum/worst case aero/"+ filename, 'var21', 'var22')
% end
% filename = "GEO_nadir_long_64";
% save("TTC_files/"+ filename)
%     end
% end
%% Plot orbit
% figure;
% plot3(mission.Satellite.TimeseriesPosICRF.Data(:,1),mission.Satellite.TimeseriesPosICRF.Data(:,2),mission.Satellite.TimeseriesPosICRF.Data(:,3))
% xlabel('x')
% ylabel('y')
% zlabel('z')
% hold on
% plot3(0, 0, 0, '*')
% moon_x = mission.SimOutput.yout{6}.Values.Data;
% %  plot3(moon_x(:,1),moon_x(:,2),moon_x(:,3))
% sun_pos = mission.SimOutput.yout{8}.Values.Data*1e-4;
% plot3([0,sun_pos(1, 1)], [0,sun_pos(1, 2)], [0,sun_pos(1, 3)])
% axis equal
% title(filename);
%savefig(['3-body-perturbation/', filename])



%% save mission
path = '/Users/mariahoflund/Library/CloudStorage/OneDrive-VinterstellarAB/Dokument/MATLAB/3-body-perturbation/';
filename = 'SSTO_no_3body';
save([path,filename]);
