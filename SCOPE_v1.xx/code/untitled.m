fversion                            =   @fluspect_bcar;

leafbio2                                        =   leafbio;

leafbio                                         =   leafbio2;
[leafopt_ref]                                   =   fversion(spectral,leafbio,optipar);
[leafopt_ref,soil]                              =   thermalspectralsignatures(options,spectral,leafbio,soil,rsfile,leafopt_ref);

% hold on; plot(spectral.wlS,leafopt.refl,'g'), xlim([400,700])
vars                                            =   {'Cab', 'Cca','Cdm','Cw','Cs','N'};
S                                               =   zeros(length(spectral.wlS),length(vars))*NaN;
for ivar = 1:length(vars)
    var                                         =   vars{ivar};
    leafbio                                     =   leafbio2;
    leafbio.(var)                               =   leafbio2.(var)*1.01;
    [leafopt]                                   =   fversion(spectral,leafbio,optipar);
    [leafopt,soil]                              =   thermalspectralsignatures(options,spectral,leafbio,soil,rsfile,leafopt);

    S(:,ivar)                                   =   abs(leafopt_ref.refl-leafopt.refl);    
end
leafbio                                         =   leafbio2;
close all

figure
semilogy(spectral.wlP,optipar.Kab,'g',spectral.wlP,optipar.Kca,'r',spectral.wlP,optipar.Kdm,'k',spectral.wlP,optipar.Kw,'b',spectral.wlP,optipar.Ks,'m'); legend('Kab','Kca','Ks','Kw','Kdm'),
xlim([400 700])

figure
plot(spectral.wlS, S(:,1),'g', spectral.wlS, S(:,2),'r--', spectral.wlS, S(:,3),'k:', spectral.wlS, S(:,4),'b.-', spectral.wlS, S(:,5),'m', spectral.wlS, S(:,6),'y>'), 
xlim([400,700])
    
legend(vars)
% legend('Cab')

leafbio                                         =   leafbio2;

