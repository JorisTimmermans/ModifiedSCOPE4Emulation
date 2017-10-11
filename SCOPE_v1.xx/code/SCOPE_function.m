% close all
function []=SCOPE_function(m)


%% Input
% input variables from ECMWF
meteo.z                 =   z;                      % [1x1]             % measurement height
meteo.Ta                 =  Ta;                     % [1x1]             % air temperature
meteo.p                 =   p;                      % [1x1]             % pressure
meteo.ea                =   RH*1;                   % [1x1]             % water vapor content
meteo.u                 =   u;                      % [1x1]             % windspeed
meteo.Rin               =   Rin;                    % [1x1]             % incoming optical radiation (not spectral)
meteo.Rli               =   Rli;                    % [1x1]             % incoming thermal radiation (not spectral)

% leaf parameters
leafbio.Cab             =   Cab;                    % [1x1]
leafbio.Cca             =   Cca;                    % [1x1]
leafbio.Cdm             =   Cdm;                    % [1x1] 
leafbio.Cw              =   Cw;                     % [1x1]
leafbio.Cs              =   Cs;                     % [1x1]
leafbio.N               =   N;                      % [1x1]
leafbio.Type            =   'C3';                   %'C3'           required for biochemical

% Canopy parameters 
canopy.LAI              =   LAI;                    % [1x1]
canopy.hc               =   hc;                     % [1x1]             % required for resistances/ stability functions
canopy.leafwidth        =   hot.*hc;                % [1x1]      %required for 1) hotspot effect (w/h) and 2) leaf boundary resistance (rbc)

%%
% solutions to hc-leafwidth 
% using Simard et al, 2011, for hc high vegetation and hc = hcmin + (hcmax-hcmin) * (NDVI - NDVImin)/(NDVImax-NDVImin), for low vegeation
% with hcmax,hcmin from % http://ldas.gsfc.nasa.gov/nldas/web/web.veg.table.htm

%% using to estimate downwelling radation in RTMo but might be simplified
atmo.M                  =   M;                      % [2162x6 double]   %atmospheric reflection/transmission/emittance 
atmo.Ta                 =   Ta;                     % [1x1]

%% set as input variables (To be provide as metadata)
angles.tts                  =   Angles(01);         % [1x1]
angles.tto                  =   Angles(02);         % [1x1]
angles.psi                  =   Angles(03);         % [1x1]

%% required by set Constant or predefined
% predefined
options.calc_ebal           =   Options(01);        % [1x1]
options.calc_vert_profiles  =   Options(02);        % [1x1]
options.calc_fluor          =   Options(03);        % [1x1]
options.calc_planck         =   Options(04);        % [1x1]
options.calc_directional    =   Options(05);        % [1x1]
options.rt_thermal          =   Options(06);        % [1x1]
options.calc_zo             =   Options(07);        % [1x1]
options.soil_heat_method    =   Options(08);        % [1x1]
options.Fluorescence_model  =   Options(09);        % [1x1]
options.calc_rss_rbs        =   Options(10);        % [1x1]
options.apply_T_corr        =   Options(11);        % [1x1]
options.verify              =   Options(12);        % [1x1]
options.save_headers        =   Options(13);        % [1x1]
options.makeplots           =   Options(14);        % [1x1]
options.simulation          =   Options(15);        % [1x1]
options.Cca_function_of_Cab =   Options(16);        % [1x1]

%%
% predefined
iter.maxit                  =   300;                % [1x1]
iter.maxEBer                =   1;                  % [1x1]
iter.Wc                     =   1;                  % [1x1]
iter.counter                =   0;                  % [1x1]

canopy.Cd                   =   0.2; 0.3; %canopy.Cd*1;        % [1x1]             %leaf drag coefficient (might be different per species?)

% Following set is constant!
canopy.nlayers              =   60;                 % [1x1]
canopy.xl                   =   linspace(0,-1,61)';                      % [61x1 double]
canopy.x                    =   canopy.xl(2:end);                       % [60x1 double]

% all interconnected, but is assumed spherical anyways!
canopy.nlazi;                   % [1x1]
canopy.lazitab;                 % [1x36 double]
canopy.nlincl;                  % [1x1]
canopy.LIDFa;                   % [1x1]
canopy.LIDFb;                   % [1x1]
canopy.hot;                     % [1x1]

% constant or defined beforehand
spectral.wlS;                   % [1x2162 double]
spectral.wlP;                   % [1x2001 double]
spectral.wlE;                   % [1x351 double]
spectral.wlF;                   % [1x211 double]
spectral.wlO;                   % [1x2001 double]
spectral.wlT;                   % [1x161 double]
spectral.wlPAR;                 % [1x301 double]
spectral.IwlP;                  % [1x2001 double]
spectral.IwlT;                  % [1x161 double]
spectral.IwlF;                  % [1x211 double]
spectral.SCOPEspec;             % [1x1 struct]
spectral.SCOPEspec.nreg;        % [1x1]
spectral.SCOPEspec.start;       % [1x3]
spectral.SCOPEspec.end;         % [1x3]
spectral.SCOPEspec.res;         % [1x3]

% required for resistances
soil.spectrum;                  % [1x1]             % should this not be defined as an option?
canopy.rb;                      % [1x1]             % leaf boundary resistance
canopy.rwc;                     % [1x1]             % within canopy layer resistance
soil.rss;                       % [1x1]              %soil resistance for evaporation (from the pore space)
soil.rbs;                       % [1x1]             %soil boundary layer resistance

%% required but might be replaced (in far future) by process-module
% soil reflectance
% soil.refl;                      % [2162x1 double]
leafbio.rho_thermal;            % 0.0100
leafbio.tau_thermal;            % 0.0100
soil.rs_thermal;                % [1x1]              %broadband soil reflectance in the thermal range

%% not Used as all?
leafbio.fqe;                    % [0.0020 0.0100]

%% Not required in Synergy framework
% not necessary
xyt.startDOY                =   NaN;                    % [1x1] 
xyt.endDOY                  =   NaN;                    % [1x1] 
xyt.LAT                     =   NaN;                    % [1x1] 
xyt.LON                     =   NaN;                    % [1x1] 
xyt.timezn                  =   NaN;                    % [1x1] 
xyt.t                       =   NaN;                    % [1x1] 
xyt.year                    =   NaN;                    % [1x1] 
xyt.k                       =   1;                      % 

% Not necessary (required to calculate soil heat flux (using force-restore), but using different G0 function)
soil.CSSOIL                 =   NaN;                    % [1x1]
soil.SMC                    =   NaN;                    % [1x1]
soil.rhos                   =   NaN;                    % [1x1]     %was 
soil.lambdas                =   NaN;                    % [1x1]
soil.cs                     =   NaN;                    % [1x1]
soil.GAM                    =   NaN;                    % [1x1]
soil.Ts                     =   [meteo.Ta,meteo.Ta]';   % [1x1]

% Not required (only for fluorescence)
% leafopt.kChlrel           =   leafopt.kChlrel*NaN;                % [2001x1 double]
% leafopt.MbI               =   leafopt.MbI*NaN;    % [211x351 double]
% leafopt.MbII              =   leafopt.MbII*NaN;  % [211x351 double]
% leafopt.MfI               =   leafopt.MfI*NaN;    % [211x351 double]
% leafopt.MfII              =   leafopt.MfII*NaN;  % [211x351 double]


canopy.CR                   =   canopy.CR*NaN;          % [1x1]             %
canopy.CD1                  =   canopy.CD1*NaN;         % [1x1]             %verhoef fitting parameter
canopy.Psicor               =   canopy.Psicor*NaN;      % [1x1]             % roughness layer correction

% required in original Biochemical.m, but not in simplified version. 
meteo.Ca                    =   NaN;                % [1x1]
meteo.Oa                    =   NaN;                % [1x1]
leafbio.Vcmo                =   NaN;                % 30
leafbio.m                   =   NaN;                % 8
canopy.kV                   =   NaN;                % [1x1]
leafbio.Rdparam             =   NaN;                % 0.0150            % dark respiration

%only used by Berry vdTol model
leafbio.Tparam              =   NaN*ones(5,1);      % [5x1 double]
        
%only used by v.Caemmerer-Magnani model
leafbio.beta                =   NaN;                % 0.5070            %
leafbio.Tyear               =   NaN;                % 15
leafbio.kNPQs               =   NaN;                % 0
leafbio.qLs                 =   NaN;                % 1
% 
leafbio.stressfactor        =   NaN;                %1   

%% Derivatives

if options.calc_zo==1
    [canopy.zo,canopy.d ]                   =   zo_and_d(soil,canopy);
elseif options.calc_zo==2
    [canopy.zo,canopy.d]                    =   zo_and_d_JT2016(canopy);
end

[canopy.lidf,canopy.litab,canopy.nlincl]    =   leafangles(canopy.LIDFa,canopy.LIDFb);

[leafopt]                                   =   fversion(spectral,leafbio,optipar);

[leafopt,soil]                              =   thermalspectralsignatures(options,spectral,leafbio,soil,rsfile,leafopt);

[rad,gap,profiles]                          =   RTMo(spectral,atmo,soil,leafopt,canopy,angles,meteo,rad,options);

%% Energy Balance Iteration
            
[iter,fluxes,rad,thermal,profiles,soil]     =   ebal(iter,options,spectral,rad,gap,                       ...
                                                leafopt,angles,meteo,soil,canopy,leafbio,xyt,k,profiles);

%                                             
% if options.calc_fluor
%     if options.calc_vert_profiles
%         [rad,profiles] = RTMf(spectral,rad,soil,leafopt,canopy,gap,angles,profiles);
%     else
%         [rad] = RTMf(spectral,rad,soil,leafopt,canopy,gap,angles,profiles);
%     end
% end

%% Post processing steps
if options.calc_planck
    rad         = RTMt_planck(spectral,rad,soil,leafopt,canopy,gap,angles,thermal.Tcu,thermal.Tch,thermal.Ts(2),thermal.Ts(1),1);
end

% if options.calc_directional
%     directional = calc_brdf(options,directional,spectral,angles,rad,atmo,soil,leafopt,canopy,meteo,profiles,thermal,leafbio);
% end

Tcu_mean                                     =  mean(mean(thermal.Tcu,1),2);

%% Plot
subplot(2,1,1,'nextplot','add')
plot(thermal.Tch,canopy.x,'Color',[0 0.5 0]) ;
plot(Tcu_mean(:),canopy.x,'Color',[0 1 0]) ;

plot(thermal.Ts(1),-1,'o','Color',[0.5 0 0]) ;
plot(thermal.Ts(2),-1,'o','Color',[1 0 0]) ;

plot(meteo.Ta,0,'o','Color',[0 0 1]) ;

legend('Shaded leaves','Sunlit leaves','Shaded soil', 'Sunlit soil','Air')  
subplot(2,1,2,'nextplot','add')
plot(spectral.wlS, rad.Lot_) ;
title(sprintf('%04.3f W/m2',rad.Lot))

%% output
Lout    =   rad.Lout;                       % [1x1]
Loutt   =   rad.Loutt;                      % [1x1]
Eoutte  =   rad.Eoutte;                     % [1x1]
LoF_    =   rad.LoF_;                       % [1x1]
Fhem_   =   rad.Fhem_;                      % [1x1]
Eouto   =   rad.Eouto;                      % [1x1]
Eout    =   rad.Eout;                       % [1x1]
Lout_   =   rad.Lout_;                      % [1x1]
Lo_     =   rad.Lo_;                        % [1x1]
PAR     =   rad.PAR;                        % [1x1]
