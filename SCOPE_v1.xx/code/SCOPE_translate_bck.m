%% Input
% input variables from ECMWF

meteo.z             =   z;                        % [1x1]             % measurement height
meteo.zm            =   zm;                       % [1x1]             % measurement height for windspeed
meteo.zh            =   zh;                       % [1x1]             % measurement height for temperature
meteo.Ta            =   Ta;                       % [1x1]             % air temperature
meteo.p             =   p;                        % [1x1]             % pressure
meteo.ea            =   ea;                       % [1x1]             % water vapor content
meteo.u             =   u;                        % [1x1]             % windspeed
meteo.Rin           =   Rin;                      % [1x1]             % incoming optical radiation (not spectral)
meteo.Rli           =   Rli;                      % [1x1]             % incoming thermal radiation (not spectral)

% leaf parameters
leafbio.Cab         =   Cab;                    % [1x1]
leafbio.Cca         =   Cca;                    % [1x1]
leafbio.Cdm         =   Cdm;                    % [1x1] 
leafbio.Cw          =   Cw;                     % [1x1]
leafbio.Cs          =   Cs;                     % [1x1]
leafbio.N           =   N;                      % [1x1]
leafbio.Vcmo        =   Vcmo;                     % [1x1]

% Canopy parameters 
canopy.LAI          =   lai;                     % [1x1]
canopy.hc           =   hc;                      % [1x1]             % required for resistances/ stability functions
canopy.leafwidth    =   lw;               % [1x1]      %required for 1) hotspot effect (w/h) and 2) leaf boundary resistance (rbc)
% canopy.leafwidth        =   hot.*hc;                % [1x1]      %required for 1) hotspot effect (w/h) and 2) leaf boundary resistance (rbc)


% using to estimate downwelling radation in RTMo but might be simplified
atmo.M(:,1)         =   M1;                         % [2162x6 double]   %atmospheric reflection/transmission/emittance 
atmo.M(:,2)         =   M2;                         % [2162x6 double]   %atmospheric reflection/transmission/emittance 
atmo.M(:,3)         =   M3;                         % [2162x6 double]   %atmospheric reflection/transmission/emittance 
atmo.M(:,4)         =   M4;                         % [2162x6 double]   %atmospheric reflection/transmission/emittance 
atmo.M(:,5)         =   M5;                         % [2162x6 double]   %atmospheric reflection/transmission/emittance 
atmo.M(:,6)         =   M6;                         % [2162x6 double]   %atmospheric reflection/transmission/emittance 

% angles
angles.tts          =   tts;
angles.tto          =   tto;
angles.psi          =   psi;

% emissivities
rho_thermal_old     =   leafbio.rho_thermal;
tau_thermal_old     =   leafbio.tau_thermal;
emissivity_leaf_old =   1 - rho_thermal_old - tau_thermal_old;

relchange           =   (1-emissivity_leaf)/(1-emissivity_leaf_old);
leafbio.rho_thermal =   rho_thermal_old * relchange;
leafbio.tau_thermal =   tau_thermal_old * relchange;

soil.rs_thermal     =   1-emissivity_soil;

clear rho_thermal_old tau_thermal_old emissivity_thermal_old


