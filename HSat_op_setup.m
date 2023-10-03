%% Start
addpath(genpath('./env_models'))
addpath(genpath('./flight_toolbox'))

%vbs_setup
%setup files
%earth_properties

% add path to
% env_models
% rvcoe,
% SimulinkToolbox,
% timetools,
% SC_orbit_dynamics,


%% Dynamics setup
%Dynamics_and_Environment/Global
ext_ENV_h0 = 700000; %600km
ext_ENV_Re = 6.37813700*10^6; %earth radius, m
ext_ENV_a0 = 3.5/log10(3); %density exponent, -
ext_ENV_rho0 = 1e-13; %air density at 600km, kg/m^3
ext_ENV_rho0 = 0; %air density at 36000km, kg/m^3
ext_ENV_solarFlux = 4.4e-6; % Solar Momentum Flux near Earth, kg/(m*s^2)
ext_ENV_muEarth = 3.986005e14; %earth gravitational constant, m^3/s^2
ext_ENV_massRatio = 0.01230002000000; %mass_moon/mass_earth
ext_ENV_massRatio_SunEarth = 3.329147856061395e+005; %mass_sun/(mass_moon+mass_earth)
ext_ENV_sunRatio = 4.652405812997204e-003; %Ratio of Radius of the Photosphere to 1 Astronomical Unit
ext_ENV_Rm = 1738200; %moon radius, m
ext_ENV_eclipseBoundary = 0.01; %boundary between penumbra and umbra, -
ext_ENV_Cs = 1358; %solar constant, W/m^2
ext_ENV_c = 2.997925E8; %light speed, m/s
ext_ENV_k_alpha_G0 = 360.9856469*pi/180; %average earth rotation, rad/day
ext_ENV_alpha_G0 = 98.8279*pi/180; %Initial Right Ascension of Greenwich Meridian, rad
ext_ENV_alpha_G0_JD = ymd2jd(1979,12,31,'gregorian')-1/2; %JD of initial RA of Greenwich Meridian
ext_ENV_phi = 109.3*pi/180; %East Longitude, rad
ext_ENV_theta = 168.6*pi/180; %Coelevation, rad
ext_ENV_Dst = 7.96826e15; %total dipole strength, Wbm
ext_ENV_AU = 149597870700; %Astronomical unit, m
ext_ENV_G= 6.674e-11; %Newtons constant ov gravitation m3⋅kg−1⋅s−2
ext_ENV_sunMass = 1.98847e30; %Mass of the sun kg
ext_ENV_moonMass = 7.342e22; %Mass of the moon kg
ext_ENV_moonDistance= 384400e3; %Distance to moon m
ext_ENV_geoRadius =    42164e3; %Radius of geostationary orbit m



if ~exist("externalsettings") | ~externalsettings
    open('HSat_orbit.slx')
end

%%EOR Setup
sim_ts_save=60;
YMD=20251201;
hms=000000;
[YYYY,MM,DD]=splitymd(YMD);
[hh,mm,ss]=splithms(hms);
sim_start = ymd2jd(YYYY,MM,DD,'gregorian')-.5+hh/24+mm/1440+ss/86400; 
[gmst]=jd2gmst(sim_start);
[sunv] = cmpsunv(YMD,hms);
rap    = atan2(sunv(2),sunv(1))+(-.5-3.5/24+hh/24+mm/1440+ss/86400)/1*2*pi;

%November 24, 4713 B.C. at noon is JD= '0
ext_Main_initDate_IC_T = sim_start;

%% Parameters for orbit raising
if ~exist("externalsettings") | ~externalsettings

    a1     = 42164e3; %m
    a1=((86164/2/pi)^2*ext_ENV_muEarth)^(1/3);
    altitude1   = a1-ext_ENV_Re;
    ecc1   = 0.0;
    p1     = a1*(1-ecc1^2); %derived
    i1     = 0*pi/180;
    ap1    = 178*pi/180;
    an1    = 180*pi/180;
%    an1    = rap-ap1;
    tran1  = 0*pi/180;
    sineccan1 = ((sqrt(1-ecc1^2)*sin(tran1))/(1+ecc1*cos(tran1)));
    eccan1 = asin(sineccan1);
    ma1    = eccan1-ecc1*sineccan1;
    ext_SC_target_COE_ECI                 = [a1,ecc1,i1,an1,ap1,tran1,eccan1,ma1]; %[sma ecc incl an ap tran];
    ext_SC_target_lambda=67*pi/180;
end
ext_SC_cutofferror_COE_ECI = [150e3, .01, 0.1*pi/180, .1*pi/180,.1*pi/180, 0, 0, 0];
ext_SC_activecorr_COE_ECI = [1, 1, 1, 0, 0, 0, 0, 0 ];
ext_SC_activecorr_COE_ECI = [1, 1, 1, 0, 0, 0, 0, 0 ];
ext_SC_balanced_COE_ECI = [ 0.10 2.10 1.40    1.0000    1.0000];

ext_EOR_atmosphere_end=1000e3; %[km]
ext_EOR_perigeeraise_end=600e3; %[km]


%sim_Target_attitude=[1     0     0     0];
%sim_Main_attitude=[0.1279   -0.8253    0.5095    0.2071];

gravity_model=3;
use_moon=1;
use_sun=1;


aero_model=0;
ext_SC_dragCoeff        = 2;
ext_SC_Normal           = [1 0 0];
ext_SC_SurfaceArea      = 1;
ext_SC_MomentArm        = [1 0 0];

use_solar_press                       = 1;
ext_SC_maxBodyArea                  = 0.7*0.8; %m^2
ext_SC_bodyReflCoeff                = 1.6;
ext_SC_SAarea                       = 0.7*0.3;
ext_SC_SAreflCoeff                  = 1.2;

%%

% enable 3-body perturbation
third_body_sun = 1;
third_body_moon = 1;
third_body_jupiter = 1;


gravity_model=3;
use_moon=1;
use_sun=1;


aero_model=0;
ext_SC_dragCoeff        = 2;
ext_SC_Normal           = [1 0 0];
ext_SC_SurfaceArea      = 1;
ext_SC_MomentArm        = [1 0 0];

use_solar_press                       = 1;
ext_SC_maxBodyArea                  = 0.7*0.8; %m^2
ext_SC_bodyReflCoeff                = 1.6;
ext_SC_SAarea                       = 0.7*0.3;
ext_SC_SAreflCoeff                  = 1.2;
%% Launch option

if ~exist("externalsettings") | ~externalsettings

    launchoption='F9GTO12.0';
    launchoption='F9SSTO12.0';
    launchoption='A6GTO6.0';

    launchoption='F9SSTO18.2';
    launchoption='F9GTO28.5';
    launchoption='F9GTO18.3';
    launchoption='A6GTO6.0';
    launchoption='COE';
    launchoption='rv';

    launchoption='A6GTO6.0';
    launchoption='COE-SFL-1';
    launchoption='F9GTO18.3';
    launchoption='F9SSTO18.2';
    launchoption='A6GTO6.0';

    launchoption='F9GTO12.0';
    launchoption='F9SSTO12.0';
    launchoption='A5GTO6.0';
    launchoption='A5SSTO6.0';

    launchoption='F9GTO28.5';
    launchoption='F9SSTO28.5';
    launchoption='F9GTO28.5';
    launchoption='A6SSTO6.0';
    launchoption='F9SSTO28.5';
    launchoption='F9SSTO17.6';
    launchoption='F9SSTO20.55';
    launchoption='rv';
    launchoption='A6GTOPG250Inc60IDV1500';
    launchoption='F9SSTO28.5';
    launchoption='A6GTOPG250inc60idv1500';
    launchoption='F9SSTO28.5';
    launchoption='F9GTOPG185inc285idv1667WCEclipse';
    launchoption='A6GTOPG250inc60idv1500WCEclipse';
    %launchoption='H3_GTOPG250inc295AG80000';
    %launchoption='GEO_67W';
    
    

    launchoption
end

switch launchoption
    case 'F9GTO12.0'
        %Ariane 6 GTO W/C alt apogee 35778 km alt perigee 256 km incl 6 deg
        %[1.2 0.9 1.1]
        alt_apogee0  = 42164e3-ext_ENV_Re;
        alt_perigee0 =  6678e3-ext_ENV_Re;
        incl0  = 12;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = 6*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;
    case 'F9SSTO12.0'
        %Falcon SSTO alt apogee 70000 km alt perigee 300 km incl 18,5 deg
        %Best sim TTO: 163 days dv 2440 m/s using 2x HEMPT: 2* 44mN
        %Best sim TTO: 189 days dv 2305 m/s using 1x PPS1350: 1* 70mN
        %[1.2 0.9 1.1]
        alt_apogee0  = 76378e3-ext_ENV_Re;
        alt_perigee0 =  6678e3-ext_ENV_Re;
        incl0  = 12;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;
        ecc0   = (ra0-rp0)/(ra0+rp0);
        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;
    case 'F9SSTO18.2'
        %Falcon SSTO alt apogee 70000 km alt perigee 300 km incl 18,5 deg
        %Best sim TTO: 163 days dv 2440 m/s using 2x HEMPT: 2* 44mN
        %Best sim TTO: 189 days dv 2305 m/s using 1x PPS1350: 1* 70mN
        alt_apogee0  = 76378e3-ext_ENV_Re;
        alt_perigee0 =  6678e3-ext_ENV_Re;
        incl0  = 18.2;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;
        ecc0   = (ra0-rp0)/(ra0+rp0);
        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;
    case 'F9SSTO28.5 old'
        %Falcon SSTO alt apogee 70000 km alt perigee 300 km incl 18,5 deg
        %Best sim TTO: 163 days dv 2440 m/s using 2x HEMPT: 2* 44mN
        %Best sim TTO: 189 days dv 2305 m/s using 1x PPS1350: 1* 70mN
        alt_apogee0  = 76378e3-ext_ENV_Re;
        alt_perigee0 =  6678e3-ext_ENV_Re;
        incl0  = 28.5;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;
        ecc0   = (ra0-rp0)/(ra0+rp0);
        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;

        ext_SC_balanced_COE_ECI = [.1 .5 .9 1 1];
    case 'F9SSTO28.5'
        %Falcon SSTO alt apogee 70000 km alt perigee 300 km incl 18,5 deg
        %Best sim TTO: 163 days dv 2440 m/s using 2x HEMPT: 2* 44mN
        %Best sim TTO: 189 days dv 2305 m/s using 1x PPS1350: 1* 70mN
        alt_apogee0  = 60000e3;
        alt_perigee0 =  185e3;
        incl0  = 28.5;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;
        ecc0   = (ra0-rp0)/(ra0+rp0);
        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;

        ext_SC_balanced_COE_ECI = [.1 .5 .9 1 1];
    case 'F9SSTO17.6'
        %Türksat 5A Wikipedia
        alt_apogee0  = 55000e3;
        alt_perigee0 =  280e3;
        incl0  = 17.6;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;
        ecc0   = (ra0-rp0)/(ra0+rp0);
        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;

        ext_SC_balanced_COE_ECI = [.1 .5 .9 1 1];
    case 'F9SSTO20.55'
        %Türksat 5A Wikipedia
        alt_apogee0  = 79341e3;
        alt_perigee0 =  185e3;
        incl0  = 20.55;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;
        ecc0   = (ra0-rp0)/(ra0+rp0);
        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;

        ext_SC_balanced_COE_ECI = [.1 .5 .9 1 1];

    case 'F9GTO28.5'
        %Falcon GTO W/C alt apogee 35786 km alt perigee 185 km incl 28,5 deg
        %Best sim TTO: 199 days dv 3016 m/s using 2x HEMPT: 2* 44mN
        %and Isp = 2200s
        alt_apogee0  = 35786000;
        alt_perigee0 =  185000;
        a0     = (alt_apogee0+alt_perigee0+2*ext_ENV_Re)/2; %24363637; %km
        ecc0=(alt_apogee0-alt_perigee0)/2/a0;
        %        ecc0   = 0.728;
        ra0    = a0*(1+ecc0);
        rp0    = a0*(1-ecc0);
        %        alt_apogee0  = ra0-ext_ENV_Re;
        %        alt_perigee0 =  rp0-ext_ENV_Re;

        p0     = a0*(1-ecc0^2); %derived
        i0     = 28.5*pi/180;
        an0    = 180*pi/180;
        ap0    = 0*pi/180;
    case 'F9GTO18.3'
        %Falcon GTO W/C alt apogee 35778 km alt perigee 256 km incl 28,5 deg'
        %Best sim TTO: almost there after 200 days dv 2451 m/s using 1x PPS1350: 70mN
        %and Isp = 2200s
        %DV used: 2450.6537 m/s
        %Did  reach orbit in 200 days
        a0     = 24396e3; %km
        ecc0   = 0.728;
        p0     = a0*(1-ecc0^2); %derived
        i0     = 18.3*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;

    case 'A5SSTO6.0'
        %Ariane 6 GTO W/C alt apogee 35778 km alt perigee 256 km incl 6 deg
        alt_apogee0  = 76378e3-ext_ENV_Re;
        alt_perigee0 =  635778628e3-ext_ENV_Re;
        incl0  = 6;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = 6*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;
    case 'A5GTO6.0'
        %Ariane 6 GTO W/C alt apogee 35778 km alt perigee 256 km incl 6 deg
        alt_apogee0  = 42164e3-ext_ENV_Re;
        alt_perigee0 =  6628e3-ext_ENV_Re;
        incl0  = 6;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = 6*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;
    case 'A6GTO6.0'
        %Ariane 6 GTO W/C alt apogee 35778 km alt perigee 256 km incl 6 deg
        alt_apogee0  = 42164e3-ext_ENV_Re;
        alt_perigee0 =  6628e3-ext_ENV_Re;
        incl0  = 6;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = 6*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;

    case 'COE-SFL-1'
        %launch orbit defined by SFL in email to Mike on June 12, 2020
        %DV used: 2098.7135 m/s
        %Reached orbit in 207 days
        alt_apogee0  = 36000e3;
        alt_perigee0 = 800e3;
        incl0  = 0;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0     = (ra0+rp0)/2;
        ecc0   = (ra0-rp0)/(ra0+rp0);
        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;

    case 'COE'
        %Ariane 6 GTO W/C alt apogee 35778 km alt perigee 256 km incl 6 deg

        a0     = 42164e3-1000e3; %km
        ecc0   = 0.728;
        ecc0   = 0.1;
        p0     = a0*(1-ecc0^2); %derived
        i0     = 1*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;

    case 'A6SSTO6.0'
        %Falcon SSTO alt apogee 70000 km alt perigee 300 km incl 18,5 deg
        %Best sim TTO: 163 days dv 2440 m/s using 2x HEMPT: 2* 44mN
        %Best sim TTO: 189 days dv 2305 m/s using 1x PPS1350: 1* 70mN
        alt_apogee0  = 60000e3;
        alt_perigee0 =  185e3;
        incl0  = 6;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;
        ecc0   = (ra0-rp0)/(ra0+rp0);
        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;

        ext_SC_balanced_COE_ECI = [.1 .5 .9 1 1];
    case 'rv'
        r                 = [8739.26e3 -41252.4e3 -5.51965e3];
        v            = [3.0078e3 0.636953e3 -0.00313559];
        [p0,a0,ecc0,i0,an0,ap0,tran0,ma0,arglat,truelon,lonper ] = rv2coe (r,v);
        %i0=7.4/2*pi/180;
        %i0=6*pi/180;
        i0=7.4*pi/180;
        an0=0;
        [r,v] = coe2rv(p0,ecc0+.002,i0,an0,ap0,tran0);
        [p0,a0,ecc0,i0,an0,ap0,tran0,ma0,arglat,truelon,lonper ] = rv2coe (r,v);
        alt_apogee0  = a0*(1+ecc0)-ext_ENV_Re;
        alt_perigee0 =  a0*(1-ecc0)-ext_ENV_Re;

    case 'A6GTOPG250inc60idv1500'
        %Ariane 6 GTO W/C 1500 m/s idv to GEO (alt apogee 35249 km) alt perigee 256 km incl 6 deg
        alt_apogee0  = 35249e3;
        alt_perigee0 =  250e3;
        incl0  = 6;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;
        if ~exist("externalsettings") | ~externalsettings

            ext_SC_balanced_COE_ECI =[ 1.0000    0.5500    0.2000    1.0000    1.0000];
        end
    case 'A6GTOPG250inc60idv1500WCEclipse'
        %Ariane 6 GTO W/C 1500 m/s idv to GEO (alt apogee 35249 km) alt perigee 256 km incl 6 deg
        alt_apogee0  = 35249e3;
        alt_perigee0 =  250e3;
        incl0  = 6;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        ap0    = 178*pi/180;
        an0    = rap-ap0;
        if ~exist("externalsettings") | ~externalsettings

            ext_SC_balanced_COE_ECI =[ 1.0000    0.5500    0.2000    1.0000    1.0000];
        end
    case 'F9GTOPG185inc285idv1800'
        %Falcon 9 GTO W/C 1700 m/s idv to GEO (alt apogee 53814 km)  alt perigee 256 km incl 28.5 deg
        alt_apogee0  = 40000e3;
        alt_perigee0 =  185e3;
        incl0  = 28.5;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;

    case 'F9GTOPG185inc285idv1700'
        %Falcon 9 GTO W/C 1700 m/s idv to GEO (alt apogee 53814 km)  alt perigee 256 km incl 28.5 deg
        alt_apogee0  = 53814e3;
        alt_perigee0 =  185e3;
        incl0  = 28.5;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;
    case 'F9GTOPG185inc285idv1667'
        %Falcon 9 GTO W/C 1700 m/s idv to GEO (alt apogee 53814 km)  alt perigee 256 km incl 28.5 deg
        alt_apogee0  = 60000e3;
        alt_perigee0 =  285e3;
        incl0  = 28.5;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;
    case 'F9GTOPG185inc285idv1600'
        %Falcon 9 GTO W/C 1700 m/s idv to GEO (alt apogee 53814 km)  alt perigee 256 km incl 28.5 deg
        alt_apogee0  = 77000e3;
        alt_perigee0 =  185e3;
        incl0  = 28.5;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;
    case 'F9SSTOPG285inc285idv1600'
        %Falcon 9 GTO W/C 1700 m/s idv to GEO (alt apogee 53814 km)  alt perigee 256 km incl 28.5 deg
        alt_apogee0  = 77000e3;
        alt_perigee0 =  185e3;
        incl0  = 28.5;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;

        ext_SC_balanced_COE_ECI = [ .500    0.70    1.100 1 1];

    case 'F9GTOPG185inc285idv1840'
        %Falcon 9 GTO W/C 1700 m/s idv to GEO (alt apogee 53814 km)  alt perigee 256 km incl 28.5 deg
        alt_apogee0  = 35700e3;
        alt_perigee0 =  185e3;
        incl0  = 28.5;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;

    case 'F9GTOPG185inc176idv1700'
        %Falcon 9 GTO W/C 1700 m/s idv to GEO (alt apogee  53000 km) alt perigee 185 km incl 17.6 deg
        alt_apogee0  = 31740e3;
        alt_perigee0 =  185e3;
        incl0  = 17.6;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;
    case 'F9GTOPG185inc176idv1500'
        %Falcon 9 GTO W/C 1500 m/s idv to GEO (alt apogee  75000 km) alt perigee 185 km incl 17.6 deg
        alt_apogee0  = 75000e3;
        alt_perigee0 =  185e3;
        incl0  = 17.6;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;

    case 'NGGTOPG250inc270idv1800'

        %NewGlenn GTO W/C 1800 m/s idv to GEO (alt apogee  25785 km) alt perigee 250 km incl 27 deg
        alt_apogee0  = 35785e3;
        alt_perigee0 =  250e3;
        incl0  = 27;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;
    case  'F9GTOPG300inc285idv1620'
        %F9 GTO W/C 1620 m/s idv to GEO (alt apogee  70000 km) alt perigee 300 km incl 28.5 deg
        alt_apogee0  = 70000e3;
        alt_perigee0 =  300e3;
        incl0  = 28.5;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;
    case '22L_SSTOPG185x80000inc285'
        %Falcon 9 GTO W/C 1700 m/s idv to GEO (alt apogee 53814 km)  alt perigee 256 km incl 28.5 deg
        alt_perigee0 =  185e3;
        alt_apogee0  = 80000e3;
        incl0  = 28.5;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;
    case 'H3_LCGTOPG2700inc201AG35586'
        %Falcon 9 GTO W/C 1700 m/s idv to GEO (alt apogee 53814 km)  alt perigee 256 km incl 28.5 deg
        alt_perigee0 =  2700e3;
        alt_apogee0  = 35586e3;
        incl0  = 20.1;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;
    case 'F9GTOPG185inc285idv1667WCEclipse'
        %Falcon 9 GTO W/C 1700 m/s idv to GEO (alt apogee 53814 km)  alt perigee 256 km incl 28.5 deg
        alt_apogee0  = 60000e3;
        alt_perigee0 =  185e3;
        incl0  = 28.5;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        ap0    = 178*pi/180;
        an0    = rap-ap0;
        an0    = 180*pi/180;
    case 'F9SSTOPG185inc285idv1600'
        %Falcon 9 GTO W/C 1700 m/s idv to GEO (alt apogee 53814 km)  alt perigee 256 km incl 28.5 deg
        alt_apogee0  = 77000e3;
        alt_perigee0 =  185e3;
        incl0  = 28.5;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;
      case 'H3_GTOPG185inc285AG35786'
        %Falcon 9 GTO W/C 1700 m/s idv to GEO (alt apogee 53814 km)  alt perigee 256 km incl 28.5 deg
        alt_perigee0 =  185e3;
        alt_apogee0  = 35786e3;
        incl0  = 28.5;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;
    H3_GTOPG185inc285AG35586      
    case 'H3_GTOPG185inc285AG35586'
        %Falcon 9 GTO W/C 1700 m/s idv to GEO (alt apogee 53814 km)  alt perigee 256 km incl 28.5 deg
        alt_perigee0 =  185e3;
        alt_apogee0  = 35586e3;
        incl0  = 28.5;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;
    
      case 'H3_SSTOPG185inc285AG80000'
        %Falcon 9 GTO W/C 1700 m/s idv to GEO (alt apogee 53814 km)  alt perigee 256 km incl 28.5 deg
        alt_perigee0 =  185e3;
        alt_apogee0  = 80000e3;
        incl0  = 28.5;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;
        
        case 'H3_GTOPG250inc150AG80000'
        %Falcon 9 GTO W/C 1700 m/s idv to GEO (alt apogee 53814 km)  alt perigee 256 km incl 28.5 deg
        alt_perigee0 =  250e3;
        alt_apogee0  = 80000e3;
        incl0  = 15.0;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 180*pi/180;
        ap0    = 180*pi/180;

case 'H3_GTOPG250inc295AG80000'
        %Falcon 9 GTO W/C 1700 m/s idv to GEO (alt apogee 53814 km)  alt perigee 256 km incl 28.5 deg
        alt_perigee0 =  250e3;
        alt_apogee0  = 80000e3;
        incl0  = 29.5;
        ra0    = alt_apogee0+ext_ENV_Re;
        rp0    = alt_perigee0+ext_ENV_Re;
        a0=(ra0+rp0)/2;

        %a0     = 24396e3; %km
        %ecc0   = 0.728;
        ecc0   = (ra0-rp0)/(ra0+rp0);

        p0     = a0*(1-ecc0^2); %derived
        i0     = incl0*pi/180;
        an0    = 0*180*pi/180;
        ap0    = 180*pi/180;

    case 'GEO_67W';
                %Falcon 9 GTO W/C 1700 m/s idv to GEO (alt apogee 53814 km)  alt perigee 256 km incl 28.5 deg
        % alt_perigee0 =  35786e3;
        % alt_apogee0  = 35786e3;
        % incl0  = 0;
        % ra0    = alt_apogee0+ext_ENV_Re;
        % rp0    = alt_perigee0+ext_ENV_Re;
        % a0=(ra0+rp0)/2;
        % 
        % %a0     = 24396e3; %km
        % %ecc0   = 0.728;
        % ecc0   = (ra0-rp0)/(ra0+rp0);
        % 
        % p0     = a0*(1-ecc0^2); %derived
        % i0     = incl0*pi/180;
        % an0    = 180*pi/180;
        % ap0    = 180*pi/180;
        a0=a1;
        ecc0=ecc1;
        i0=i1;
        an0=an1;
        ap0=ap1;
        tran0=tran1;
        eccan0=eccan1;
        ma0=ma1;
        p0     = a0*(1-ecc0^2); %derived
        otherwise
        disp('Launch orbit not found. Discontinuing')
        return
end
%% 
tran0  = 0*pi/180;

sineccan0 = ((sqrt(1-ecc0^2)*sin(tran0))/(1+ecc0*cos(tran0)));
eccan0 = asin(sineccan0);
ma0    = eccan0-ecc0*sineccan0;

arglat = ap0+tran0;
lonper =an0+ap0;
truelon=lonper+tran0;

ext_SC_init_COE_ECI                 = [a0,ecc0,i0,an0,ap0,tran0,eccan0,ma0]; %[sma ecc incl an ap tran];
ext_SC_init_COE_ECI(ext_SC_init_COE_ECI-ext_SC_target_COE_ECI==0)=1e-16;
[r,v] = coe2rv(p0,ecc0,i0,an0,ap0,tran0,arglat,truelon,lonper);
r=r;
v=v;
ext_SC_init_pos_ECI                 = [r];
ext_SC_init_vel_ECI                 = [v];
%ext_SC_init_pos_ECI                 = [8739.26e3 -41252.4e3 -5.51965e3];
%ext_SC_init_momentum_ECI            = [3.0078e3 0.636953e3 -0.00313559]*sim_SC_mass;

%%

%% Spacecraft setup
if ~exist("externalsettings") | ~externalsettings

    %spacecraft mass
    SC_conf='espaRingLimited';
    %SC_conf='SSL-1300';
    SC_conf='OR40mN';
    SC_conf='MMIOR';
    SC_conf='SSL-1300';
    SC_conf='Program A';
    SC_conf='OR80mN';
    SC_conf='OR800mN';
    SC_conf='OR400mN';
    SC_conf='OR200mN';
    SC_conf='OR80mN';
    SC_conf='OR90mN';
    SC_conf='OR110mN';
    %SC_conf='OR550mN';
    SC_conf='NewBaseline';
end

switch SC_conf
    case "espaRingLimited"
        sim_SC_init_wet_mass=535; %Kg
        sim_SC_init_propellant_mass=300; %Kg
        ORthrusternumber=2;
        ORthrustermodel='HEMPT';
    case 'SSL-1300'
        sim_SC_init_wet_mass=5000; %Kg
        sim_SC_init_propellant_mass=300; %Kg
        ORthrusternumber=4;
        ORthrustermodel='SPT-140';

    case 'OR400mN'
        sim_SC_init_wet_mass=560; %Kg
        sim_SC_init_propellant_mass=300; %Kg
        ORthrusternumber=10;
        ORthrustermodel='SPT-70';
    case 'OR200mN'
        sim_SC_init_wet_mass=560; %Kg
        sim_SC_init_propellant_mass=300; %Kg
        ORthrusternumber=5;
        ORthrustermodel='SPT-70';
    case 'OR80mN'
        sim_SC_init_propellant_mass=110; %Kg
        sim_SC_init_wet_mass=355+52+sim_SC_init_propellant_mass; %Kg
        ORthrusternumber=2;
        ORthrustermodel='SPT-70';
    case 'OR800mN'
        sim_SC_init_wet_mass=560; %Kg
        sim_SC_init_propellant_mass=300; %Kg
        ORthrusternumber=20;
        ORthrustermodel='SPT-70';
    case 'OR40mN'
        sim_SC_init_wet_mass=350; %Kg
        sim_SC_init_propellant_mass=300; %Kg
        ORthrusternumber=40;
        ORthrustermodel='1mN';
    case 'OR90mN'
        sim_SC_init_wet_mass=400; %Kg
        sim_SC_init_propellant_mass=300; %Kg
        ORthrusternumber=2;
        ORthrustermodel='SPT-70';
    case 'OR110mN'
        sim_SC_init_wet_mass=700; %Kg
        sim_SC_init_propellant_mass=300; %Kg
        ORthrusternumber=2;
        ORthrustermodel='OSE ORTh';
    case 'OR550mN'
        sim_SC_init_wet_mass=700; %Kg
        sim_SC_init_propellant_mass=300; %Kg
        ORthrusternumber=2;
        ORthrustermodel='OSE ORThx5';
    case 'OR25mN'
        sim_SC_init_wet_mass=350; %Kg
        sim_SC_init_propellant_mass=300; %Kg
        ORthrusternumber=25;
        ORthrustermodel='1mN';
    case 'MMIOR'
        sim_SC_init_wet_mass=560; %Kg
        sim_SC_init_propellant_mass=300; %Kg
        ORthrusternumber=4;
        ORthrustermodel='MMI';
    case 'Program A'
        sim_SC_init_wet_mass=560; %Kg
        sim_SC_init_propellant_mass=300; %Kg
        ORthrusternumber=1;
        ORthrustermodel='PPS1350';

    case 'NewBaseline'
        if ~exist("externalsettings") | ~externalsettings
            sim_SC_init_wet_mass=934; %Kg
        end
        sim_SC_init_propellant_mass=207; %Kg
        sim_SC_init_dry_mass=sim_SC_init_wet_mass-sim_SC_init_propellant_mass;
        ORthrusternumber=2;

        ORthrustermodel='BHT600(900W)';
        SKthrustermodel='EV0';


    case 'NewBaseline_1xEOR'
        if ~exist("externalsettings") | ~externalsettings
            sim_SC_init_wet_mass=750; %Kg
        end
        sim_SC_init_propellant_mass=207; %Kg
        ORthrusternumber=1;
        OR_cant_angle=0*pi/180;
        ORthrustermodel='BHT600(900W)';
        SKthrustermodel='EV0';

    otherwise
        disp('Configuration not found. Discontinuing')
        return


end


sim_SC_dry_mass=sim_SC_init_wet_mass-sim_SC_init_propellant_mass; %Kg

%% RTU and attitude guidance setup
%Reduced mass
ext_ENV_muEarth = 3.986005e14; %earth gravitational constant, m^3/s^2
%J2 coefficient
ext_ENV_J2=0.00108263000000;
%Earth radiuse
GNC_glb_Re_s = 6.37813700*10^6; %earth radius, m




%Time step [seconds]
ts=1;


%Thrust profile
Motor_active=[
    0     0 0 0 0];

%Target
q_slo2scb_tgt=[1 0 0 0];



switch SKthrustermodel
    case 'EV0'
        SK_thrusterforce=13e-3;
        %Thrust effective exhaust velocity [m/s; Ns/kg] ==G0*Isp
        %        OR_cant_angle=17.5*pi/180;
        SK_Isp=1335;
        SK_veff=SK_Isp*9.81; %SPT-70
    case 'BHT350'
end

switch ORthrustermodel
    case 'HEMPT'

        %Thrust force [N]

        OR_force=ORthrusternumber*44e-3;

        %Thrust effective exhaust velocity [m/s; Ns/kg] ==G0*Isp

        OR_veff=22000; %PPS1350

        %Thrust mass flow [kg/s]
    case 'PPS1350'
        OR_force=ORthrusternumber*90e-3;
        %Thrust effective exhaust velocity [m/s; Ns/kg] ==G0*Isp

        OR_veff=16600; %PPS1350

    case 'SPT-140'
        OR_force=ORthrusternumber*300e-3;
        %Thrust effective exhaust velocity [m/s; Ns/kg] ==G0*Isp

        OR_veff=15000; %PPS1350

    case 'SPT-70'
        OR_force=ORthrusternumber*45e-3;
        %Thrust effective exhaust velocity [m/s; Ns/kg] ==G0*Isp

        OR_veff=1840*9.81; %SPT-70

    case 'OSE ORTh'
        OR_thrusterforce=55e-3;
        OR_force=ORthrusternumber*OR_thrusterforce;
        %Thrust effective exhaust velocity [m/s; Ns/kg] ==G0*Isp
        OR_cant_angle=17.5*pi/180;
        OR_Isp=1840;
        OR_veff=1840*9.81; %SPT-70
        OR_veffGrav=OR_veff*cos(OR_cant_angle); %SPT-70 %After gravity losses
    case 'BHT600(900W)'
        OR_thrusterforce=53e-3;

        %Thrust effective exhaust velocity [m/s; Ns/kg] ==G0*Isp
        %        OR_cant_angle=17.5*pi/180;
        OR_Isp=1597;
        OR_veff=OR_Isp*9.81; %SPT-70
    case 'OSE ORThx5'
        OR_thrusterforce=5*55e-3;
        OR_force=ORthrusternumber*OR_thrusterforce;
        %Thrust effective exhaust velocity [m/s; Ns/kg] ==G0*Isp
        OR_cant_angle=17.5*pi/180;
        OR_Isp=1840;
        OR_veff=1840*9.81; %SPT-70
        OR_veffGrav=OR_veff*cos(OR_cant_angle); %SPT-70 %After gravity losses



    case '1mN'
        OR_force=ORthrusternumber*1e-3;
        %Thrust effective exhaust velocity [m/s; Ns/kg] ==G0*Isp

        OR_veff=40000; %PPS1350

    case 'MMI'
        OR_force=ORthrusternumber*20e-3;
        %Thrust effective exhaust velocity [m/s; Ns/kg] ==G0*Isp

        OR_veff=20000; %MMI


    otherwise
        disp('Launch orbit not found. Discontinuing')
        return


end
OR_force=ORthrusternumber*OR_thrusterforce;
OR_active=[0 0 0  0 0 0 0 0 ORthrusternumber>=1 ORthrusternumber>=3 ORthrusternumber>=2];
%Thrust mass flow [kg/s]


%Direction of acceleration

%from 20221224
sim_SC_init_MoI_EOR_wet=[
    1789.4664 -7.6657 14.0713
    -7.6657 815.5687 -2.1297
    14.0713 -2.1297 1862.6698];

OR_vec=[0 0  1];
sim_init_propellant=207;
sim_EOR_propellant=139.29;
sim_SC_init_CoG_propellant=[0 0 0.24758];
sim_SC_init_CoG_EOR_dry=[-0.00448 0.00106 0.46931];
sim_SC_init_CoG_EOR_wet=(sim_SC_init_CoG_EOR_dry*sim_SC_init_dry_mass+sim_SC_init_CoG_propellant*sim_SC_init_propellant_mass)/sim_SC_init_wet_mass;
sim_SC_mid_CoG_EOR_wet=(sim_SC_init_CoG_EOR_dry*sim_SC_init_dry_mass+(sim_init_propellant-sim_EOR_propellant/2)*sim_SC_init_CoG_propellant)/(sim_SC_init_wet_mass+(sim_init_propellant-sim_EOR_propellant/2));

sim_ESK_position=[
    -0.65  -0.5 -0.1+sim_SC_init_CoG_EOR_dry(3)
    -0.65  -0.5 +0.1+sim_SC_init_CoG_EOR_dry(3)
    +0.65  -0.5 -0.1+sim_SC_init_CoG_EOR_dry(3)
    +0.65  -0.5 +0.1+sim_SC_init_CoG_EOR_dry(3)
    +0.65  +0.5 -0.1+sim_SC_init_CoG_EOR_dry(3)
    +0.65  +0.5 +0.1+sim_SC_init_CoG_EOR_dry(3)
    -0.65  +0.5 -0.1+sim_SC_init_CoG_EOR_dry(3)
    -0.65  +0.5 +0.1+sim_SC_init_CoG_EOR_dry(3)
    ];
sim_ESK_thrustvec_nom=[
    +1  +1  0
    +1  +1  0
    -1  +1  0
    -1  +1  0
    -1  -1  0
    -1  -1  0
    +1  -1  0
    +1  -1  0
    ]/sqrt(2);

sim_ESK_thrustvec_nom_thrusterf=[
    0 0 1
    0 0 1
    0 0 1
    0 0 1
    0 0 1
    0 0 1
    0 0 1
    0 0 1];
sim_ESK_scb2thruster_nom_eul=[
    0 +90 +45
    0 +90 +45
    0 -90 -45
    0 -90 -45
    0 -90 +45
    0 -90 +45
    0 +90 -45
    0 +90 -45
    ]*pi/180;
sim_ESK_scb2thruster_nom_quat=eul2quat(fliplr(sim_ESK_scb2thruster_nom_eul));


sim_ESK_thrustvec_err_eul=[
    0 1 0
    0 1 0
    0 1 0
    0 1 0
    0 1 0
    0 1 0
    0 1 0
    0 1 0
    ]*pi/180*0;
sim_ESK_thrustvec_err_quat=eul2quat(fliplr(sim_ESK_thrustvec_err_eul));

%sim_ESK_thrustvec_nom=quatrotate(quatinv(sim_ESK_scb2thruster_nom_quat),sim_ESK_thrustvec_nom_thrusterf)
%sim_ESK_thrustvec_n=quatrotate(quatinv(sim_ESK_thrustvec_err_quat),sim_ESK_thrustvec_nom);

sim_ESK_thrustvec_n_thrusterf=quatrotate(quatinv(sim_ESK_thrustvec_err_quat),sim_ESK_thrustvec_nom_thrusterf);
sim_ESK_thrustvec_n=quatrotate(quatinv(sim_ESK_scb2thruster_nom_quat),sim_ESK_thrustvec_n_thrusterf);


%OR_cant_angle=17.5*pi/180;
%        EOR_zposition=-.050;
EOR_zposition=-.180; %Number from CAD on Nov 23, 2022
%  EOR_rposition=tan(30*pi/180)*.075+.050;
EOR_rposition=0.1525;%Number from CAD on Nov 23, 2022

OR_cant_angle=atan(-EOR_rposition/(EOR_zposition-sim_SC_init_CoG_EOR_wet(3)));

sim_EOR_position=[
    EOR_rposition 0 EOR_zposition
    [cos(2*pi/3) sin(2*pi/3)]*EOR_rposition  EOR_zposition
    [cos(4*pi/3) sin(4*pi/3)]*EOR_rposition  EOR_zposition
    ];


%         sim_EOR_thrustvec_nom=[
%             -sin(OR_cant_angle) 0 cos(OR_cant_angle)
%             - [cos(2*pi/3) sin(2*pi/3)]*sin(OR_cant_angle) cos(OR_cant_angle)
%             - [cos(4*pi/3) sin(4*pi/3)]*sin(OR_cant_angle) cos(OR_cant_angle)
%             ];

sim_EOR_thrustvec_nom_thrusterf=[
    0 0 1
    0 0 1
    0 0 1];
sim_EOR_scb2thruster_nom_eul=[
    0 -OR_cant_angle*180/pi 0
    0 -OR_cant_angle*180/pi 120
    0 -OR_cant_angle*180/pi 240
    ]*pi/180;
sim_EOR_scb2thruster_nom_quat=eul2quat(fliplr(sim_EOR_scb2thruster_nom_eul));

%         sim_EOR_thrustvec_err_eul=[
%             1 1 0
%             1 1 0
%             1 1 0
%             ]*pi/180*0;
sim_EOR_thrustvec_err_eul=[
    1.0 -0 0
    1.0 -0 0
    1.0 -0 0
    ]*pi/180*0;

sim_EOR_thrustvec_err_eul=sim_EOR_thrustvec_err_eul;

sim_EOR_thrustvec_err_quat=eul2quat(fliplr(sim_EOR_thrustvec_err_eul));

%sim_EOR_thrustvec_err_quat=eul2quat(fliplr(sim_EOR_thrustvec_err_eul));
%sim_EOR_thrustvec_n=quatrotate(quatinv(sim_EOR_thrustvec_err_quat),sim_EOR_thrustvec_nom);

%         sim_EOR_thrustvec_nom=quatrotate(quatinv(sim_EOR_scb2thruster_nom_quat),sim_EOR_thrustvec_nom_thrusterf);
%         sim_EOR_thrustvec_n=quatrotate(quatinv(sim_EOR_thrustvec_err_quat),sim_EOR_thrustvec_n_thrusterf);

sim_EOR_thrustvec_n_thrusterf=quatrotate(quatinv(sim_EOR_thrustvec_err_quat),sim_EOR_thrustvec_nom_thrusterf);
sim_EOR_thrustvec_n=quatrotate(quatinv(sim_EOR_scb2thruster_nom_quat),sim_EOR_thrustvec_n_thrusterf);


sim_thruster_position=[
    sim_ESK_position
    sim_EOR_position];

sim_thrust_vector=[
    sim_ESK_thrustvec_n
    sim_EOR_thrustvec_n];

sim_thrust_force=[
    SK_thrusterforce
    SK_thrusterforce
    SK_thrusterforce
    SK_thrusterforce
    SK_thrusterforce
    SK_thrusterforce
    SK_thrusterforce
    SK_thrusterforce
    OR_thrusterforce
    OR_thrusterforce
    OR_thrusterforce
    ];

sim_thrust_massflow=[
    SK_thrusterforce/SK_veff
    SK_thrusterforce/SK_veff
    SK_thrusterforce/SK_veff
    SK_thrusterforce/SK_veff
    SK_thrusterforce/SK_veff
    SK_thrusterforce/SK_veff
    SK_thrusterforce/SK_veff
    SK_thrusterforce/SK_veff
    OR_thrusterforce/OR_veff
    OR_thrusterforce/OR_veff
    OR_thrusterforce/OR_veff
    ];



