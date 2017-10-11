clc
close all
clear all
fclose all;

%% Initialize SCOPE

SCOPE
leafbio.Cdm=0.0001;

NwlTt                           =   length(rad.Lot_);

%%
% meteo.Ta  = 10;
% meteo.Rin = 900;
% meteo.Rli = 400;
SCOPE_run

figure('Position',[20 20 1024 800])
plot(thermal.Ts(1),-1,'rd',thermal.Ts(2),-1,'or',thermal.Tch,canopy.x,'>g',permute(mean(mean(thermal.Tcu,1),2),[3 1 2]),canopy.x,'*g',meteo.Ta,0.1,'bo'); 

% hold on
% meteo.Rin = 100;
% SCOPE_run

% plot(thermal.Ts(1),-1,'rd',thermal.Ts(2),-1,'or',thermal.Tch,canopy.x,'>g',permute(mean(mean(thermal.Tcu,1),2),[3 1 2]),canopy.x,'*g',meteo.Ta,0.1,'bo'); 
legend('Tsh','Tsu','Tch','Tsu','Ta')
title('Temperature profile for different components')


Min                                 =   [ 0.1, 0.01,  0.1, 0.001, 0.01, 1.0, 0.01, 0.01];
Max                                 =   [60.0, 0.15, 40.0, 0.040, 8.00, 4.0, 0.10, 0.05];

%%
varnames                            =   {'Cab','Cw','Cca','Cdm', 'Cs',  'N', 'rho_thermal','tau_thermal'};
Nvar                                =   length(varnames);

LAI                                 =   logspace(log10(0.01),log10(5),10);
Rin                                 =   linspace(50,900,6);

% individual run for different prospect parameters
S_Lot_                              =      zeros(NwlTt,Nvar,length(LAI), length(Rin),5)*NaN;
S_Ts                                =   zeros(2,Nvar,length(LAI), length(Rin),5)*NaN;
[S_Tch,S_Tcu]                       =   deal(zeros(canopy.nlayers,Nvar,length(LAI), length(Rin),5)*NaN);


for i = 1:length(Rin)
    meteo.Rin                       =   Rin(i);
%     legendstr{iRin}                 =   sprintf('Rin = % 3.1fW/m2',Rin(iRin));
    for ilai = 1:length(LAI)
        % LAI sensitivity reference run 1
        canopy.LAI                  =   LAI(ilai)*1.01;
        SCOPE_run
        Lot_ref_                    =   rad.Lot_;
        Ts_ref                      =   thermal.Ts;
        Tch_ref                     =   thermal.Tch;
        Tcu_ref                     =   permute(mean(mean(thermal.Tcu,1),2),[3 1 2]);
        
        % LAI sensitivity reference run 1
        canopy.LAI                  =   LAI(ilai);
        SCOPE_run
        S_Lot_ref(:,ilai)           =   abs(rad.Lot_ - Lot_ref_);
        S_Ts_ref(:,ilai)            =   abs(thermal.Ts - Ts_ref);
        S_Tch_ref(:,ilai)           =   abs(thermal.Tch - Tch_ref);
        S_Tcu_ref(:,ilai)           =   abs(permute(mean(mean(thermal.Tcu,1),2),[3 1 2]) - Tcu_ref);
        
        
        % reference run
        leafbio_ref                 =   leafbio;
        canopy_ref                  =   canopy;
        %%
        for j = 1:Nvar
            varname                 =   varnames{j};
            Val                     =   linspace(Min(j),Max(j),5);
            for jj=1:length(Val)
                % reset variables
                leafbio             =   leafbio_ref;
                leafbio.(varname)   =   Val(jj);
                

                SCOPE_run
                
                Lot_                =   rad.Lot_;
                Ts                  =   thermal.Ts;
                Tch                 =   thermal.Tch;
                Tcu                 =   permute(mean(mean(thermal.Tcu,1),2),[3 1 2]);
                
                % infi
                leafbio.(varname)   =   Val(jj)*1.01;
                SCOPE_run       
                
                S_Lot_s(:,jj)       =   abs(rad.Lot_ - Lot_);
                S_Tss(:,jj)         =   abs(thermal.Ts - Ts);
                S_Tchs(:,jj)        =   abs(thermal.Tch - Tch);
                S_Tcus(:,jj)        =   abs(permute(mean(mean(thermal.Tcu,1),2),[3 1 2]) - Tcu);
            end
            
            leafbio             =   leafbio_ref;                

    %             Ts_ref = 0;
    %             Tcu_ref = 0;
    %             Tch_ref = 0;
%                 S_Lot_(:,j,ilai,i)   =   abs(Lot_ref_ - rad.Lot_);
%                 S_Ts(:,j,ilai,i)     =   abs(Ts_ref   - thermal.Ts);
%                 S_Tch(:,j,ilai,i)    =   abs(Tch_ref  - thermal.Tch);
%                 S_Tcu(:,j,ilai,i)    =   abs(Tcu_ref  - permute(mean(mean(thermal.Tcu,1),2),[3 1 2]));

                S_Lot_(:,j,ilai,i,:)   =   S_Lot_s;
                S_Ts(:,j,ilai,i,:)     =   S_Tss;
                S_Tch(:,j,ilai,i,:)    =   S_Tchs;
                S_Tcu(:,j,ilai,i,:)    =   S_Tcus;
        end
    end
end


%% Post Processing
% MODIS bands
[~,iband_1]                         =   min((spectral.wlS - (10.780 + 11.280)/2*1e3).^2);
[~,iband_2]                         =   min((spectral.wlS - (11.770 + 12.270)/2*1e3).^2);

S_Lot_1                             =   permute(mean(S_Lot_(iband_1,:,:,:,:),5),[3 4 2 1]);
S_Lot_2                             =   permute(mean(S_Lot_(iband_2,:,:,:,:),5),[3 4 2 1]);     %sunlit

S_Tsh                               =   permute(mean(S_Ts(1,:,:,:,:),5),[3 4 2 1]);
S_Tsu                               =   permute(mean(S_Ts(2,:,:,:,:),5),[3 4 2 1]);     %sunlit


S_Tch_t                             =   permute(mean(S_Tch( 1,:,:,:,:),5),[3 4 2 1]);
S_Tcu_t                             =   permute(mean(S_Tcu( 1,:,:,:,:),5),[3 4 2 1]);   %sunlit


S_Tch_m                             =   permute(mean(S_Tch(30,:,:,:,:),5),[3 4 2 1]);
S_Tcu_m                             =   permute(mean(S_Tcu(30,:,:,:,:),5),[3 4 2 1]);   %sunlit


S_Tch_b                             =   permute(mean(S_Tch(60,:,:,:,:),5),[3 4 2 1]);
S_Tcu_b                             =   permute(mean(S_Tcu(60,:,:,:,:),5),[3 4 2 1]);   %sunlit

%% Quantiles
q                                   =   [0.05, 0.5, 0.95];
for iq=1:length(q)
    S_Lot_1b(:,iq,:)                 =   quantile(S_Lot_1,q(iq),2);
    S_Lot_2b(:,iq,:)                 =   quantile(S_Lot_2,q(iq),2);

    S_Tshb(:,iq,:)                   =   quantile(S_Tsh,q(iq),2);
    S_Tsub(:,iq,:)                   =   quantile(S_Tsu,q(iq),2);


    S_Tch_tb(:,iq,:)                 =   quantile(S_Tch_t,q(iq),2);
    S_Tcu_tb(:,iq,:)                 =   quantile(S_Tcu_t,q(iq),2);


    S_Tch_mb(:,iq,:)                 =   quantile(S_Tch_m,q(iq),2);
    S_Tcu_mb(:,iq,:)                 =   quantile(S_Tcu_m,q(iq),2);


    S_Tch_bb(:,iq,:)                 =   quantile(S_Tch_b,q(iq),2);
    S_Tcu_bb(:,iq,:)                 =   quantile(S_Tcu_b,q(iq),2);
end
legendstr = {' 5%%','50%%','95%%'};


S_Lot_1                 =   S_Lot_1b;
S_Lot_2                 =   S_Lot_2b;

S_Tsh                   =   S_Tshb;
S_Tsu                   =   S_Tsub;


S_Tch_t                 =   S_Tch_tb;
S_Tcu_t                 =   S_Tcu_tb;


S_Tch_m                 =   S_Tch_mb;
S_Tcu_m                 =   S_Tcu_mb;


S_Tch_b                 =   S_Tch_bb;
S_Tcu_b                 =   S_Tcu_bb;

%%
Nparam = length(varnames);
Nc = 3;
Nr = ceil(Nparam/Nc);

% max(S_Tch_b,[],3)
%% Soil
figure('Position',[20 20 1024 800])
for j=1:Nparam
    h(j) = subplot(Nr,Nc,j);
%     semilogy(LAI,max(S_Tsu,[],3),'color',[0.8 0.8 0.8]); hold on
    semilogy(LAI,S_Tsu(:,:,j)+1e-12);
    ylabel(varnames{j})
    
    ylim([1e-6 1e0])
    set(gca,'ytick',[logspace(-6,1,6)])
    grid on    
    
    if j==2
        title('Prospect sensitivity of Sunlit Soil Temperature (for different LAI values)')
    end
    if j==5
        legend(legendstr)
    end
    if j>6
        xlabel('LAI [m2/m2]')
    end
end



figure('Position',[20 20 1024 800])
for j=1:Nparam
    h(j) = subplot(Nr,Nc,j);
    semilogy(LAI,S_Tsh(:,:,j)+1e-12); hold on
    ylabel(varnames{j})
    
    
    ylim([1e-5 1e0])
    set(gca,'ytick',[logspace(-5,1,6)])
    grid on    
    
    if j==2
        title('Prospect sensitivity of Shaded Soil Temperature (for different LAI values)')
    end
    
    if j==5
        legend(legendstr)
    end
    if j>6
        xlabel('LAI [m2/m2]')
    end
end


%% Vegetation
figure('Position',[20 20 1024 800])
for j=1:Nparam
    h(j) = subplot(Nr,Nc,j);
    semilogy(LAI,S_Tcu_b(:,:,j));
    ylabel(varnames{j})
    
    ylim([1e-5 1e0])
    set(gca,'ytick',[logspace(-5,1,6)])
    grid on    
    
    if j==2
        title('Prospect sensitivity of Sunlit Bottom Canopy Temperature (for different LAI values)')
    end
    if j==5
        legend(legendstr)
    end
    if j>6
        xlabel('LAI [m2/m2]')
    end
end


figure('Position',[20 20 1024 800])
for j=1:Nparam
    h(j) = subplot(Nr,Nc,j);
    semilogy(LAI,S_Tch_b(:,:,j));
    ylabel(varnames{j})
    
    ylim([1e-5 1e0])
    set(gca,'ytick',[logspace(-5,1,6)])
    grid on    
    
    if j==2
        title('Prospect sensitivity of Shaded Bottom Canopy Temperature (for different LAI values)')
    end
    if j==5
        legend(legendstr)
    end
    if j>6
        xlabel('LAI [m2/m2]')
    end
end


figure('Position',[20 20 1024 800])
for j=1:Nparam
    h(j) = subplot(Nr,Nc,j);
    semilogy(LAI,S_Tcu_t(:,:,j));
    ylabel(varnames{j})
    
    ylim([1e-5 1e0])
    set(gca,'ytick',[logspace(-5,1,6)])
    grid on    
    
    if j==2
        title('Prospect sensitivity of Sunlit Top Canopy Temperature (for different LAI values)')
    end
    if j==5
        legend(legendstr)
    end
    if j>6
        xlabel('LAI [m2/m2]')
    end
end



figure('Position',[20 20 1024 800])
for j=1:Nparam
    h(j) = subplot(Nr,Nc,j);
    semilogy(LAI,S_Tch_t(:,:,j));
    ylabel(varnames{j})
    
    ylim([1e-5 1e0])
    set(gca,'ytick',[logspace(-5,1,6)])
    grid on    
    
    if j==2
        title('Prospect sensitivity of Shaded Top Canopy Temperature (for different LAI values)')
    end
    if j==5
        legend(legendstr)
    end
    if j>6
        xlabel('LAI [m2/m2]')
    end
end


%% Radiation
figure('Position',[20 20 1024 800])
for j=1:Nparam
    h(j) = subplot(Nr,Nc,j);
    semilogy(LAI,S_Lot_1(:,:,j));
    ylabel(varnames{j})
    
    ylim([1e-6 1e0])
    set(gca,'ytick',[logspace(-6,1,6)])
    grid on    
    
    if j==2
        title('Prospect sensitivity of Outgoing Radiation @ MODIS Band31 (for different LAI values)')
    end
    if j==5
        legend(legendstr)
    end
    if j>6
        xlabel('LAI [m2/m2]')
    end
end


figure('Position',[20 20 1024 800])
for j=1:Nparam
    h(j) = subplot(Nr,Nc,j);
    semilogy(LAI,S_Lot_2(:,:,j));
    ylabel(varnames{j})
    
    ylim([1e-6 1e0])
    set(gca,'ytick',[logspace(-6,1,6)])
    grid on    
    
    if j==2
        title('Prospect sensitivity of Outgoing Radiation @ MODIS Band32 (for different LAI values)')
    end
    if j==5
        legend(legendstr)
    end
    if j>6
        xlabel('LAI [m2/m2]')
    end
end


