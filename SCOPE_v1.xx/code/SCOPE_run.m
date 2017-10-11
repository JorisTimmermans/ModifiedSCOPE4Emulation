% close all


% if exist('leafbio2','var')
%     leafbio = leafbio2;
% end
% if exist('canopy2','var')
%     canopy = canopy2;
% end
% if exist('meteo2','var')    
%     meteo = meteo2; 
% end
%% Input
% % input variables from ECMWF
% meteo.z;                        % [1x1]             % measurement height
% meteo.Ta;                       % [1x1]             % air temperature
% meteo.p;                        % [1x1]             % pressure
% meteo.ea;                       % [1x1]             % water vapor content
% meteo.u;                        % [1x1]             % windspeed
% meteo.Rin;                      % [1x1]             % incoming optical radiation (not spectral)
% meteo.Rli;                      % [1x1]             % incoming thermal radiation (not spectral)
% 
% % leaf parameters
% leafbio.Cab;                    % [1x1]
% leafbio.Cca;                    % [1x1]
% leafbio.Cdm;                    % [1x1] 
% leafbio.Cw;                     % [1x1]
% leafbio.Cs;                     % [1x1]
% leafbio.N;                      % [1x1]
% 
% % Canopy parameters 
% canopy.LAI;                     % [1x1]
% canopy.hc;                      % [1x1]             % required for resistances/ stability functions
% canopy.leafwidth;               % [1x1]      %required for 1) hotspot effect (w/h) and 2) leaf boundary resistance (rbc)
% % canopy.leafwidth        =   hot.*hc;                % [1x1]      %required for 1) hotspot effect (w/h) and 2) leaf boundary resistance (rbc)
% 
% %%
% % solutions to hc-leafwidth 
% % using Simard et al, 2011, for hc high vegetation and hc = hcmin + (hcmax-hcmin) * (NDVI - NDVImin)/(NDVImax-NDVImin), for low vegeation
% % with hcmax,hcmin from % http://ldas.gsfc.nasa.gov/nldas/web/web.veg.table.htm
% 
% %% using to estimate downwelling radation in RTMo but might be simplified
% atmo.M;                         % [2162x6 double]   %atmospheric reflection/transmission/emittance 
% atmo.Ta;                        % [1x1]
% 
% 
%        
% 
% %% required by set Constant or predefined
% % predefined
% options.calc_ebal;              % [1x1]
% options.calc_vert_profiles;     % [1x1]
% options.calc_fluor;             % [1x1]
% options.calc_planck;            % [1x1]
% options.calc_directional;       % [1x1]
% options.rt_thermal;             % [1x1]
% options.calc_zo;                % [1x1]
% options.soil_heat_method;       % [1x1]
% options.Fluorescence_model;     % [1x1]
% options.calc_rss_rbs;           % [1x1]
% options.apply_T_corr;           % [1x1]
% options.verify;                 % [1x1]
% options.save_headers;           % [1x1]
% options.makeplots;              % [1x1]
% options.simulation;             % [1x1]
% options.Cca_function_of_Cab;    % [1x1]
% 
% canopy.Cd=canopy.Cd*1;          % [1x1]             %leaf drag coefficient (might be different per species?)
% 
% % predefined
% iter.maxit;                     % [1x1]
% iter.maxEBer;                   % [1x1]
% iter.Wc;                        % [1x1]
% iter.counter;                   % [1x1]
% 
% % set as input variables (To be provide as metadata)
% angles.tts;                     % [1x1]
% angles.tto;                     % [1x1]
% angles.psi;                     % [1x1]
% 
% 
% % Following set is constant!
% canopy.nlayers;                 % [1x1]
% canopy.x;                       % [60x1 double]
% canopy.xl;                      % [61x1 double]
% 
% 
% % all interconnected, but is assumed spherical anyways!
% canopy.nlazi;                   % [1x1]
% canopy.lazitab;                 % [1x36 double]
% canopy.nlincl;                  % [1x1]
% canopy.LIDFa;                   % [1x1]
% canopy.LIDFb;                   % [1x1]
% canopy.hot;                     % [1x1]
% 
% % constant or defined beforehand
% spectral.wlS;                   % [1x2162 double]
% spectral.wlP;                   % [1x2001 double]
% spectral.wlE;                   % [1x351 double]
% spectral.wlF;                   % [1x211 double]
% spectral.wlO;                   % [1x2001 double]
% spectral.wlT;                   % [1x161 double]
% spectral.wlPAR;                 % [1x301 double]
% spectral.IwlP;                  % [1x2001 double]
% spectral.IwlT;                  % [1x161 double]
% spectral.IwlF;                  % [1x211 double]
% spectral.SCOPEspec;             % [1x1 struct]
% spectral.SCOPEspec.nreg;        % [1x1]
% spectral.SCOPEspec.start;       % [1x3]
% spectral.SCOPEspec.end;         % [1x3]
% spectral.SCOPEspec.res;         % [1x3]
% 
% % required for resistances
% soil.spectrum;                  % [1x1]             % should this not be defined as an option?
% canopy.rb;                      % [1x1]             % leaf boundary resistance
% canopy.rwc;                     % [1x1]             % within canopy layer resistance
% soil.rss;                       % [1x1]              %soil resistance for evaporation (from the pore space)
% soil.rbs;                       % [1x1]             %soil boundary layer resistance
% 
% %% required but might be replaced (in far future) by process-module
% % soil reflectance
% % soil.refl;                      % [2162x1 double]
% leafbio.rho_thermal;            % 0.0100
% leafbio.tau_thermal;            % 0.0100
% soil.rs_thermal;                % [1x1]              %broadband soil reflectance in the thermal range
% 
% %% not Used as all?
% leafbio.fqe;                    % [0.0020 0.0100]
% 
% %% Not required in Synergy framework
% % not necessary
% xyt.k                                       = 1;                      % note necessary
% xyt.startDOY                                =   xyt.startDOY*NaN;  % [1x1] 
% xyt.endDOY                                  =   xyt.endDOY*NaN;      % [1x1] 
% xyt.LAT                                     =   xyt.LAT*NaN;            % [1x1] 
% xyt.LON                                     =   xyt.LON*NaN;            % [1x1] 
% xyt.timezn                                  =   xyt.timezn*NaN;      % [1x1] 
% xyt.t                                       =   xyt.t*NaN;                % [1x1] 
% xyt.year                                    =   xyt.year*NaN;          % [1x1] 
% 
% % Not necessary (required to calculate soil heat flux (using force-restore), but using different G0 function)
% soil.CSSOIL                                 =   soil.CSSOIL*NaN;    % [1x1]
% soil.SMC                                    =   soil.SMC*NaN;        % [1x1]
% soil.rhos                                   =   soil.rhos*NaN;      % [1x1]     %was 
% soil.lambdas                                =   soil.lambdas*NaN;% [1x1]
% soil.cs                                     =   soil.cs*NaN;                        % [1x1]
% soil.GAM                                    =   soil.GAM*NaN;        % [1x1]
% soil.Ts                                     =   [meteo.Ta,meteo.Ta]'; % [1x1]
% 
% % Not required (only for fluorescence)
% % leafopt.kChlrel                           =   leafopt.kChlrel*NaN;                % [2001x1 double]
% % leafopt.MbI                               =   leafopt.MbI*NaN;    % [211x351 double]
% % leafopt.MbII                              =   leafopt.MbII*NaN;  % [211x351 double]
% % leafopt.MfI                               =   leafopt.MfI*NaN;    % [211x351 double]
% % leafopt.MfII                              =   leafopt.MfII*NaN;  % [211x351 double]
% 
% 
% % ?
% canopy.CR                                   =   canopy.CR*NaN;        % [1x1]             %
% canopy.CD1                                  =   canopy.CD1*NaN;      % [1x1]             %verhoef fitting parameter
% canopy.Psicor                               =   canopy.Psicor*NaN;% [1x1]             % roughness layer correction
% 
% % required in original Biochemical.m, but not in simplified version. 
% meteo.Ca                                    =   meteo.Ca*NaN;          % [1x1]
% meteo.Oa                                    =   meteo.Oa *NaN;        % [1x1]
% leafbio.Type;                   %'C3'           required for biochemical
% leafbio.Vcmo                                =   leafbio.Vcmo*NaN;  % 30
% leafbio.m                                   =   leafbio.m*NaN;          % 8
% canopy.kV                                   =   canopy.kV*NaN;        % [1x1]
% leafbio.Rdparam                             =   leafbio.Rdparam*NaN;                % 0.0150            % dark respiration
% 
% %only used by Berry vdTol model
% leafbio.Tparam                              =   leafbio.Tparam*NaN;                 % [5x1 double]
%         
% %only used by v.Caemmerer-Magnani model
% leafbio.beta                                =   leafbio.beta*NaN;                      % 0.5070            %
% leafbio.Tyear                               =   leafbio.Tyear*NaN;                    % 15
% leafbio.kNPQs                               =   leafbio.kNPQs*NaN;                    % 0
% leafbio.qLs                                 =   leafbio.qLs*NaN;                        % 1
% % 
% leafbio.stressfactor                        =   leafbio.stressfactor*NaN;      %1   

%%
% zh = -10;
% zm = -2;
% % 
% Ta = 40.926339181278507;
% p = 1048.7344094844721;
% ea = 87.563631668792993;
% u = 4.6133797261068006;
% Rin = 363.66906227891866;
% Rli = 477.89406381864916;
% Cab = 33.600471716199593;
% Cca = 0.75275256047436767;
% Cdm = 0.0066919337382604482;
% Cw = 0.01516847704208084;
% Cs = 0.24508933768651955;
% N = 1.0121204011473792;
% LAI = 0.020908142236247401;
% hc = 48.470652317810654;
% lw = 0.076361861620745375;
% 
% % 
% SCOPE_translate_bck
  
% meteo             =   rmfield(meteo,'z');

%%

% meteo.Ta = 35;
% canopy

if options.calc_zo==1
    [canopy.zo,canopy.d ]                   =   zo_and_d(soil,canopy);
elseif options.calc_zo==2
    [canopy.zo,canopy.d]                    =   zo_and_d_JT2016(canopy);
end

[canopy.lidf,canopy.litab,canopy.nlincl]    =   leafangles(canopy.LIDFa,canopy.LIDFb);


[leafopt]                                   =   fversion(spectral,leafbio,optipar);
[leafopt,soil]                              =   thermalspectralsignatures(options,spectral,leafbio,soil,rsfile,leafopt);

% use transmissitivities/reflectivities so that the output is consistent with MCD43B
% leafopt.refl(spectral.IwlPAR)               =   leaf_refPAR;
% leafopt.tran(spectral.IwlPAR)               =   leaf_transPAR;
% soil.refl(spectral.IwlPAR)                  =   soil_refPAR;
% 
% leafopt.refl(spectral.IwlNIR)               =   leaf_refNIR;
% leafopt.tran(spectral.IwlNIR)               =   leaf_transNIR;
% soil.refl(spectral.IwlNIR)                  =   soil_refNIR;


[rad,gap,profiles]                          =   RTMo(spectral,atmo,soil,leafopt,canopy,angles,meteo,rad,options);

[iter,fluxes,rad,thermal,profiles,soil]     =   ebal(iter,options,spectral,rad,gap,                       ...
                                                leafopt,angles,meteo,soil,canopy,leafbio,xyt,k,profiles);

% if options.calc_fluor
%     if options.calc_vert_profiles
%         [rad,profiles]                      =   RTMf(spectral,rad,soil,leafopt,canopy,gap,angles,profiles);
%     else
%         [rad]                               =   RTMf(spectral,rad,soil,leafopt,canopy,gap,angles,profiles);
%     end
% end

% profile on, 
if options.calc_planck
    rad                                     =   RTMt_planck(spectral,rad,soil,leafopt,canopy,gap,angles,thermal.Tcu,thermal.Tch,thermal.Ts(2),thermal.Ts(1),1);
end

% profile report
if options.calc_directional
    directional                             =   calc_brdf(options,directional,spectral,angles,rad,atmo,soil,leafopt,canopy,meteo,profiles,thermal,leafbio);
end

rad.Ltot_                                   =   rad.Lo_ + rad.Lot_;

%% Band simulation
[MODIS.TEB, LUT_error_min]                  =   Simulate_MODIS_TEB(MODIS, rad, spectral);


% plot(rad.Lot_)
% keyboard
% profile off
%%
% hold on
% plot(spectral.wlS, rad.rdd,...
%      spectral.wlS, rad.rsd,...
%      spectral.wlS, rad.rdo,...
%      spectral.wlS, rad.rso)
% legend('dd','sd','do','so')
% xlim([0.4, 3000])

