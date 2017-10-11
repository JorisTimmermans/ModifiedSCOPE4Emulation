%% 14. Run the model
fprintf('\n The calculations start now \r')
calculate = 1;

for k = 1:telmax
    
    if options.simulation == 1, 
        vi(vmax>1)                          =   k; 
    end
    
    if options.simulation == 0, 
        vi(vmax==telmax)                    =   k; 
    end
    
    [soil,leafbio,canopy,meteo,angles,xyt]  =   select_input(V,vi,canopy,options,xyt,soil);
    
    % added by JT to enable different measurement heights for windspeed and temperature
    try
        meteo.zm                            =   meteo.z;
        meteo.zh                            =   meteo.z;
    end
        
    % perform quality checks
    if options.simulation ~=1,
        fprintf('simulation %i ', k );
        fprintf('of %i \n', telmax);
    else        
        calculate = 0;
        if k>=I_tmin && k<=I_tmax
            quality_is_ok                   =   ~isnan(meteo.p*meteo.Ta*meteo.ea*meteo.u.*meteo.Rin.*meteo.Rli);
            fprintf('time = %4.2f \n', xyt.t(k));
            if quality_is_ok
                calculate = 1;
            end
        end
    end
    
    if calculate        
        iter.counter                        =   0;
        canopy.lidf                         =   leafangles(canopy.LIDFa,canopy.LIDFb);    % This is 'ladgen' in the original SAIL model,
        [canopy.lidf,...
         canopy.litab,...
         canopy.nlincl]                     =   leafangles(canopy.LIDFa,canopy.LIDFb);    % This is 'ladgen' in the original SAIL model,
        
        fversion                            =   @fluspect_bcar;
%         [leafopt] = fversion(spectral,leafbio,optipar);
        
%         IwlP                              =   spectral.IwlP;
%         IwlT                              =   spectral.IwlT;
%         
%         rho(IwlP)                         =   leafopt.refl;
%         tau(IwlP)                         =   leafopt.tran;
%         rlast                             =   rho(nwlP);
%         tlast                             =   tau(nwlP);
%         
%         rs(IwlP)                          =   soil.refl;
%         
%         rslast                            =   rs(nwlP);
%         
%         switch options.rt_thermal
%             case 0
%                 rho(IwlT)                 =   ones(nwlT,1) * leafbio.rho_thermal;
%                 tau(IwlT)                 =   ones(nwlT,1) * leafbio.tau_thermal;
%                 rs(IwlT)                  =   ones(nwlT,1) * soil.rs_thermal;
%             case 1
%                 rho(IwlT)                 =   ones(nwlT,1) * rlast;
%                 tau(IwlT)                 =   ones(nwlT,1) * tlast;
%                 rs(IwlT)                  =   ones(nwlT,1) * rslast;
%         end
%         leafopt.refl                      =   rho;     % extended wavelength ranges are stored in structures
%         leafopt.tran                      =   tau;
%         soil.refl                         =   rs;
        
        soil.Ts                             =   meteo.Ta * ones(2,1);       % initial soil surface temperature
        
        if length(F(4).FileName)>1 && options.simulation==0
            atmfile                         =   [path_input 'radiationdata/' char(F(4).FileName(k))];
            atmo.M                          =   aggreg(atmfile,spectral.SCOPEspec);
        end
        atmo.Ta                             =   meteo.Ta;
        
        

        if options.calc_ebal==1
            
            SCOPE_translate_fwd
            SCOPE_translate_bck            
            SCOPE_run
            return
                
        else
            [rad,gap,profiles]              =   RTMo(spectral,atmo,soil,leafopt,canopy,angles,meteo,rad,options);
            
            Fc                              =   (1-gap.Ps(1:end-1))'/nl;      %           Matrix containing values for Ps of canopy
            
            fluxes.aPAR                     =   canopy.LAI*(Fc*rad.Pnh        + meanleaf(canopy,rad.Pnu    , 'angles_and_layers',gap.Ps));% net PAR leaves
            fluxes.aPAR_Cab                 =   canopy.LAI*(Fc*rad.Pnh_Cab    + meanleaf(canopy,rad.Pnu_Cab, 'angles_and_layers',gap.Ps));% net PAR leaves
            fluxes.aPAR_Wm2                 =   canopy.LAI*(Fc*rad.Rnh_PAR    + meanleaf(canopy,rad.Rnu_PAR, 'angles_and_layers',gap.Ps));% net PAR leaves
            
            if options.calc_fluor
                profiles.etah               =   ones(60,1);
                profiles.etau               =   ones(13,36,60);
                if options.calc_vert_profiles
                    [rad,profiles]          =   RTMf(spectral,rad,soil,leafopt,canopy,gap,angles,profiles);
                else
                    [rad]                   =   RTMf(spectral,rad,soil,leafopt,canopy,gap,angles,profiles);
                end
            end
                
                
                
        end
        output_data
    end
    if options.simulation==2 && telmax>1, vi  = count(nvars,vi,vmax,1); end
end

% rad.rsd
% rad.rdd
% rad.rdo
% rad.rso
% 
% rad.Lot
% rad.Lo_
