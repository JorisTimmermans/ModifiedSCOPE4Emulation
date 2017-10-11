% varnames            =   fieldnames(rad);
% for j=1:length(varnames)
%     varname         =   varnames{j};
%     value           =   rad.(varname);
%     ierror          =   isnan(value);
%     value(ierror)   =   -9999;
%     rad.(varname)   =   value;
% end


%% output requried for retrieval
rsd     =   rad.rsd;                            % [nwl x 1 double]
rdd     =   rad.rdd;                            % [nwl x 1 double]
rdo     =   rad.rdo;                            % [nwl x 1 double]
rso     =   rad.rso;                            % [nwl x 1 double]

Lo_     =   rad.Lo_;                            % [nwl x 1 double]      %directional radiance
Lot_    =   rad.Lot_;
Ltot_   =   rad.Ltot_;

Tb      =   MODIS.TEB;
Bandnr  =   MODIS.nr;
Bandwl  =   (MODIS.min+MODIS.max)/2;

%% Output 
Ts      =   thermal.Ts;
Tcu     =   thermal.Tcu;
Tch     =   thermal.Tch;

Rntot   =   fluxes.Rntot;
lEtot   =   fluxes.lEtot;
Htot    =   fluxes.Htot;
Gtot    =   fluxes.Gtot;

% %%
% 
% Eoutte  =   rad.Eoutte;                         % 430.7214
% Eouto   =   rad.Eouto;                          % 104.2612
% 
% vb      =   rad.vb;                             % [2162x1 double]
% vf      =   rad.vb;                             % [2162x1 double]
% 
% 
% %% used as input or intermediary.. so not required as output
% Esun_   =   rad.Esun_;                          % [2162x1 double]       %not 
% Esky_   =   rad.Esky_;                          % [2162x1 double]
% 
% fEsuno  =   rad.fEsuno;                         % [2162x1 double]
% fEskyo  =   rad.fEskyo;                         % [2162x1 double]
% fEsunt  =   rad.fEsunt;                         % [2162x1 double]
% fEskyt  =   rad.fEskyt;                         % [2162x1 double]
% 
% Rnhc    =   rad.Rnhc;                           % [60x1 double]
% Rnuc    =   rad.Rnuc;                           % [13x36x60 double]
% 
% Pnh     =   rad.Pnh;                            % [60x1 double]
% Pnu     =   rad.Pnu;                            % [13x36x60 double]
% Pnh_Cab =   rad.Pnh_Cab;                        % [60x1 double]
% Pnu_Cab =   rad.Pnu_Cab;                        % [13x36x60 double]
% Rnh_PAR =   rad.Rnh_PAR;                        % [60x1 double]
% Rnu_PAR =   rad.Rnu_PAR;                        % [13x36x60 double]
% 
% %% not defined
% Lout    =   rad.Lout;                           % [1x1]
% Loutt   =   rad.Loutt;                          % [1x1]
% LoF_    =   rad.LoF_;                           % [1x1]
% Fhem_   =   rad.Fhem_;                          % [1x1]
% Eout    =   rad.Eout;                           % [1x1]
% Lout_   =   rad.Lout_;                          % [1x1]
% 
% PAR     =   rad.PAR;                            % [1x1]
% inPAR   =   rad.inPAR;                          % 0.0013
% 
% Eplu_   =   rad.Eplu_;                          % [61x2162 double]
% Emin_   =   rad.Emin_;                          % [61x2162 double]
% 
% Eout_   =   rad.Eout_;                          % [2162x1 double]
% 
% Eoutt   =   rad.Eoutt;                          % 1.5541
% Rnhs    =   rad.Rnhs;                           % 73.6959
% Rnus    =   rad.Rnus;                           % 411.3289
% Etoto   =   rad.Etoto;                          % 645.6923
% 
% Emint   =   rad.Emint;                          % [61x1 double]
% Eplut   =   rad.Eplut;                          % [61x1 double]
% 
% Lot     =   rad.Lot;                            % 140.0434
% Rnust   =   rad.Rnust;                          % -138.3805
% Rnhst   =   rad.Rnhst;                          % -23.6031
% LotBB_  =   rad.LotBB_;                         % [2162x1 double]
% Lot_    =   rad.Lot_;                           % [2162x1 double]
% Eoutte_ =   rad.Eoutte_;                        % [2162x1 double]
% Eplut_  =   rad.Eplut_;                         % [61x161 double]
% Emint_  =   rad.Emint_;                         % [61x161 double]
% 
% Rnuct   =   rad.Rnuct;                          % [13x36x60 double]
% Rnhct   =   rad.Rnhct;                          % [60x1 double]
