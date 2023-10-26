function [alt_apogee0, alt_perigee0, incl0, ra0, rp0, a0, ecc0, p0, i0, an0, ap0, ext_SC_balanced_COE_ECI] = ...
    launch_settings(launchoption, ext_ENV_Re)

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

        ext_SC_balanced_COE_ECI =[ 5.0000    2.5500    0.5000    1.0000    1.0000];
        
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

    case 'F9SSTOPG250inc285idv1600'
        %Falcon 9 GTO W/C 1700 m/s idv to GEO (alt apogee 53814 km)  alt perigee 256 km incl 28.5 deg
        alt_apogee0  = 35249e3;
        alt_perigee0 =  250e3;
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

        ext_SC_balanced_COE_ECI = [ .1000 1.9000 3.3000 1 1];

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