function [leafopt,soil] = thermalspectralsignatures(options,spectral,leafbio,soil,rsfile,leafopt)
IwlP                = spectral.IwlP;
IwlT                = spectral.IwlT;

nwlP                =   length(IwlP);
nwlT                =   length(IwlT);
[rho,tau,rs] = deal(zeros(nwlP + nwlT,1));


rho(IwlP)           = leafopt.refl;
tau(IwlP)           = leafopt.tran;


rlast               = rho(nwlP);
tlast               = tau(nwlP);


rs(IwlP)            = rsfile(:,soil.spectrum+1);
rslast              = rs(nwlP);

switch options.rt_thermal
    case 0
        rho(IwlT)   = ones(nwlT,1) * leafbio.rho_thermal;
        tau(IwlT)   = ones(nwlT,1) * leafbio.tau_thermal;
        rs(IwlT)    = ones(nwlT,1) * soil.rs_thermal;
    case 1
        rho(IwlT)   = ones(nwlT,1) * rlast;
        tau(IwlT)   = ones(nwlT,1) * tlast;
        rs(IwlT)    = ones(nwlT,1) * rslast;
end
leafopt.refl = rho;     % extended wavelength ranges are stored in structures
leafopt.tran = tau;
soil.refl    = rs;
