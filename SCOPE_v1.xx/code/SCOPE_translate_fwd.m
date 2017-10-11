%% Input
% input variables from ECMWF
z       =   meteo.z;                        % [1x1]             % measurement height
zm      =   meteo.zm;                       % [1x1]             % measurement height for windspeed
zh      =   meteo.zh;                       % [1x1]             % measurement height for temperature
Ta      =   meteo.Ta;                       % [1x1]             % air temperature
p       =   meteo.p;                        % [1x1]             % pressure
ea      =   meteo.ea;                       % [1x1]             % water vapor content
u       =   meteo.u;                        % [1x1]             % windspeed
Rin     =   meteo.Rin;                      % [1x1]             % incoming optical radiation (not spectral)
Rli     =   meteo.Rli;                      % [1x1]             % incoming thermal radiation (not spectral)

% z, Ta, p, ea, u, Rin, Rli

% leaf parameters
Cab     =   leafbio.Cab;                    % [1x1]
Cca     =   leafbio.Cca;                    % [1x1]
Cdm     =   leafbio.Cdm;                    % [1x1] 
Cw      =   leafbio.Cw;                     % [1x1]
Cs      =   leafbio.Cs;                     % [1x1]
N       =   leafbio.N;                      % [1x1]
Vcmo    =   leafbio.Vcmo;                   % [1x1]

% Cab,Cca, Cdm, Cw, Cs, N

% Canopy parameters 
lai     =   canopy.LAI;                     % [1x1]
hc      =   canopy.hc;                      % [1x1]             % required for resistances/ stability functions
lw      =   canopy.leafwidth;               % [1x1]      %required for 1) hotspot effect (w/h) and 2) leaf boundary resistance (rbc)
% canopy.leafwidth        =   hot.*hc;                % [1x1]      %required for 1) hotspot effect (w/h) and 2) leaf boundary resistance (rbc)
% LAI, hc, lw

% using to estimate downwelling radation in RTMo but might be simplified
M1      =   atmo.M(:,1);                         % [2162x6 double]   %atmospheric reflection/transmission/emittance 
M2      =   atmo.M(:,2);                         % [2162x6 double]   %atmospheric reflection/transmission/emittance 
M3      =   atmo.M(:,3);                         % [2162x6 double]   %atmospheric reflection/transmission/emittance 
M4      =   atmo.M(:,4);                         % [2162x6 double]   %atmospheric reflection/transmission/emittance 
M5      =   atmo.M(:,5);                         % [2162x6 double]   %atmospheric reflection/transmission/emittance 
M6      =   atmo.M(:,6);                         % [2162x6 double]   %atmospheric reflection/transmission/emittance 

% angles
tts     =   angles.tts;
tto     =   angles.tto;
psi     =   angles.psi;


% emissivities
emissivity_leaf     = 1 - leafbio.rho_thermal - leafbio.tau_thermal;
emissivity_soil     = 1 - soil.rs_thermal;


% rho_thermal = leafbio.rho_thermal
% tau_thermal = leafbio.tau_thermal
% soil.rs_thermal


