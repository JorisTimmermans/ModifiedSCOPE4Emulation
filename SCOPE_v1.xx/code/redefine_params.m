function [leafbio,meteo,canopy,atmo,rad] =   redefine_params(leafbio,meteo,canopy,atmo,Vars,rad,i)

global constants
constants                       =   define_constants;

%%
fprintf('\n')
varnames                        =   fieldnames(Vars)';
Nvar                            =   length(varnames);
Nsample                         =   zeros(1,Nvar);
for ivar=1:Nvar
    name                        =   varnames{ivar};
    
%     fprintf(name)
%     fprintf(' %f ',Vars.(name))
    
    Nsample(ivar)               =   length(Vars.(name));
%     fprintf(' %f\n',Nsample(ivar))
end

% we run this code simulateously for the sensitivity analysis and the temporal evolution to create priors. This satisfies that both can run parallel
% of each other
try
    Vars.Cab;
    [I.Cab,I.Cw,I.Cca,I.Cdm,I.Cs,I.N,I.rho_thermal, I.tau_thermal, I.LAI, I.Ta, I.Rin, I.Rli, I.p, I.ea, I.u]          =   ind2sub(Nsample, i);
catch
    [I.LAI,I.hc,I.zh,I.zm,I.Ta,I.ea,I.p,I.u,I.Rin,I.Rli]          =   deal(i); %ind2sub(Nsample, i);
end

% keyboard

rad_names                       =   fieldnames(rad);
leafbio_names                   =   fieldnames(leafbio);
meteo_names                     =   fieldnames(meteo);
canopy_names                    =   fieldnames(canopy);
atmo_names                      =   fieldnames(atmo);

for ivar=1:Nvar
    name                        =   varnames{ivar};
    ii                          =   I.(name);
    
    
    vars                        =   Vars.(name);
    

    
    switch name
        case rad_names
            rad.(name)          =   vars(ii);            
    end            
    switch name
        case leafbio_names
            leafbio.(name)      =   vars(ii);            
    end
    switch name
        case meteo_names
            meteo.(name)        =   vars(ii);
    end
    switch name
        case canopy_names
            canopy.(name)       =   vars(ii);                
    end
    switch name
        case atmo_names        
            atmo.(name)         =   vars(ii);
    end
       
end