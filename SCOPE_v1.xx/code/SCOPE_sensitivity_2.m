clc
close all
clear all;
fclose all;

%% Parameters
varnames                            =   {'Cab', 'Cw',   'Cca',  'Cdm',  'Cs',  'N',     'rho_thermal',  'tau_thermal','LAI',    'Ta',   'Rin',  'Rli',  'p',    'ea', 'u'};
Min                                 =   [ 0.1,  0.01,   0.1,    0.001,  0.01,   1.0,    0.01,           0.01         ,0.001,    0.00,     0.0,  100,     800.0,   8,  0.01];
Max                                 =   [60.0,  0.15,   40.0,   0.040,  8.00,   4.0,    0.10,           0.05         ,7.000,    50.0,   900.0,  500,    1100.0,  22,  8.00];
filename                            =   'LSP';
Nsample                             =   [4,        4,      4,       4,    4,      4,       4,              4,             4,       1,       3,    3,         1,   1,     1]; 

filename                            =   'Meteo4';
Nsample                             =   [1,        1,      1,       1,    1,      1,       1,              1,             4,       4,       4,    4,         4,   4,     4]; 

% filenames for restart mechanism  
outputdir                           =   '../output/Sensitivity_Analysis/';
path2file                           =   [outputdir,filename,'-sensitivity-Simulations-Intermediate.mat'];
path2file2                          =   [outputdir,filename,'-sensitivity-Simulations-Finished.mat'];





%% Initialize SCOPE

SCOPE

leafbio.Cdm                         =   0.0001;

NwlTt                               =   length(rad.Lot_);

% Prerun SCOPE
% meteo.Ta  = 10;
% meteo.Rin = 900;
% meteo.Rli = 400;
SCOPE_run

% figure('Position',[20 20 1024 800])
% plot(thermal.Ts(1),-1,'rd',thermal.Ts(2),-1,'or',thermal.Tch,canopy.x,'>g',permute(mean(mean(thermal.Tcu,1),2),[3 1 2]),canopy.x,'*g',meteo.Ta,0.1,'bo'); 

% hold on
% meteo.Rin = 100;
% SCOPE_run

% plot(thermal.Ts(1),-1,'rd',thermal.Ts(2),-1,'or',thermal.Tch,canopy.x,'>g',permute(mean(mean(thermal.Tcu,1),2),[3 1 2]),canopy.x,'*g',meteo.Ta,0.1,'bo'); 
% legend('Tsh','Tsu','Tch','Tsu','Ta')
% title('Temperature profile for different components')

%% Derived Parameters
% Nsample                             =   min(Nsample,1);
% Nsample(1:4)                        =   2;

Nvar                                =   length(varnames);
% Nsample                             =   5;
Vars                                =   struct('Cab',[], 'Cw',[], 'Cca',[],  'Cdm',[] , 'Cs',[],  'N',[],     'rho_thermal',[],  'tau_thermal',[],'LAI',[],    'Ta',[],   'Rin',[],  'Rli',[] , 'p',[],    'ea',[], 'u',[]);
for ivar = 1:Nvar
    name                            =   varnames{ivar};
    Vars.(name)                     =   linspace(Min(ivar), Max(ivar),Nsample(ivar));
%     x(ivar,:)                       =   Vars.(name);
end

Nsamples                            =   prod(Nsample);
% sizeV                               =   Nsample;

iwls                                =   2013:5:length(rad.Lot_); %plot(spectral.wlS(iwls),rad.Lot_(iwls)); 
ilayer                              =   unique([1,1:10:60,60]);

nwl                                 =   length(iwls);
nlayer                              =   length(ilayer);

% Ts                                  =   zeros([   2,Nsample],'uint16');
% Tch                                 =   zeros([  nlayer,Nsample],'uint16')*NaN;
% Tcu                                 =   zeros([  nlayer,Nsample],'uint16')*NaN;
% Lot_                                =   zeros([nwl,Nsample],'uint16')*NaN;

%% Allocate Matrices
Ts                                  =   zeros([   2,Nsamples],'uint16');
Tch                                 =   zeros([  nlayer,Nsamples],'uint16')*NaN;
Tcu                                 =   zeros([  nlayer,Nsamples],'uint16')*NaN;
Lot_                                =   zeros([nwl,Nsamples],'uint16')*NaN;

%% Derived Parameters
% Nsample                             =   min(Nsample,1);
% Nsample(1:4)                        =   2;

Nvar                                =   length(varnames);
% Nsample                             =   5;
Vars                                =   struct('Cab',[], 'Cw',[], 'Cca',[],  'Cdm',[] , 'Cs',[],  'N',[],     'rho_thermal',[],  'tau_thermal',[],'LAI',[],    'Ta',[],   'Rin',[],  'Rli',[] , 'p',[],    'ea',[], 'u',[]);
for ivar = 1:Nvar
    name                            =   varnames{ivar};
    Vars.(name)                     =   linspace(Min(ivar), Max(ivar),Nsample(ivar));
%     x(ivar,:)                       =   Vars.(name);
end

Nsamples                            =   prod(Nsample);
% sizeV                               =   Nsample;

iwls                                =   2013:5:length(rad.Lot_); %plot(spectral.wlS(iwls),rad.Lot_(iwls)); 
ilayer                              =   unique([1,1:10:60,60]);

nwl                                 =   length(iwls);
nlayer                              =   length(ilayer);

% Ts                                  =   zeros([   2,Nsample],'uint16');
% Tch                                 =   zeros([  nlayer,Nsample],'uint16')*NaN;
% Tcu                                 =   zeros([  nlayer,Nsample],'uint16')*NaN;
% Lot_                                =   zeros([nwl,Nsample],'uint16')*NaN;




%% Processing
duration                            =   [];
p                                   =   gcp();
tile                                =   25;
if exist(path2file2,'file')
    load(path2file2)
else
    if exist(path2file,'file')
        load(path2file)
        close(h)
        istart                          =   iii;
        gcp
    else
        istart                          =   1;
    end

    %%
    counter                             =   0;
    h = waitbar(0,'Please wait...');
    for iii=istart:tile:Nsamples
        tic

        ratio                           =   (iii-1)/Nsamples;
        string                          =   [sprintf('Please wait,..%3.0f-%3.0f   (%3.1f%%)',iii, Nsamples, ratio*100), duration];
        waitbar(i/Nsamples,h,string)



        %%
        t_start                            =   (iii+0);
        t_end                              =   min(iii+(tile-1),Nsamples);
%         parfor i = t_start:t_end
        parfor i = t_start:t_end    
            % Set model-parameters  
            [leafbio2,meteo2,canopy2,atmo2, rad2]          =   redefine_params(leafbio,meteo,canopy,atmo,Vars,rad,i);

            % Run the model 
            [Rad,Thermal]                                   =   SCOPE_run_as_function(leafopt,leafbio2,optipar, canopy2, soil, rsfile, atmo2, meteo2, options, spectral, angles, iter, xyt, fversion, rad2,k, MODIS);




        % store output
            Thermal.Tcu_mean        =   mean(mean(Thermal.Tcu,1),2);
            Ts(:,i)                 = uint16(Thermal.Ts*100);
            Tch(:,i)                = uint16(Thermal.Tch(ilayer)*100);
            Tcu(:,i)                = uint16(Thermal.Tcu_mean(ilayer)*100);
            Lot_(:,i)               = uint16(Rad.Lot_(iwls)*100);

        %     Lot_s(:,jj)       =   abs(rad.Lot_ - Lot_);
        %     Tss(:,jj)         =   abs(thermal.Ts - Ts);
        %     Tchs(:,jj)        =   abs(thermal.Tch - Tch);
        %     Tcus(:,jj)        =   abs(permute(mean(mean(thermal.Tcu,1),2),[3 1 2]) - Tcu);

        end

        dt                          =   toc;
        samples2compute             =   (Nsamples-iii)/tile;
        duration                    =   sprintf(',= %03.1f days remaining',dt*samples2compute / 60/60/24);

        if counter==40
            save(path2file)
            delete(gcp)
            gcp
            counter                 =   0;
        end
        counter                     =   counter+1;
    end
    waitbar(1,h,sprintf('Finished'))
    close(h)
    save(path2file2)

    
end
return

%% Post analysis
% load data
load(path2file2)

%% Perform Sensitivity Analysis
% Reshape data in appropriate form
Tsh                             =   reshape(single(Ts(1,:))/100, Nsample);  % shaded soil
Tsu                             =   reshape(single(Ts(2,:))/100, Nsample);  % sunlit soil

Tch_t                           =   reshape(single(Tch(1,:))/100, Nsample); % shaded top leaves
Tch_m                           =   reshape(single(Tch(4,:))/100, Nsample); % shaded medium leaves
Tch_b                           =   reshape(single(Tch(7,:))/100, Nsample); % shaded bottom leaves

Tcu_t                           =   reshape(single(Tcu(1,:))/100, Nsample); % sunlit top leaves
Tcu_m                           =   reshape(single(Tcu(4,:))/100, Nsample); % sunlit top leaves
Tcu_b                           =   reshape(single(Tcu(7,:))/100, Nsample); % sunlit top leaves

% sensitivity analysis
[Vname,Vmean(1,:)]                   =   SensitivityAnalysis(Vars,Tsh,filename,'Shaded Soil Temperature',1);
[Vname,Vmean(2,:)]                   =   SensitivityAnalysis(Vars,Tsu,filename,'Sunlit Soil Temperature',1);

[Vname,Vmean(3,:)]                   =   SensitivityAnalysis(Vars,Tch_t,filename,'Shaded Top Leaf Temperature',1);
[Vname,Vmean(4,:)]                   =   SensitivityAnalysis(Vars,Tch_m,filename,'Shaded Middle Leaf Temperature',1);
[Vname,Vmean(5,:)]                   =   SensitivityAnalysis(Vars,Tch_b,filename,'Shaded Bottom Leaf Temperature',1);

[Vname,Vmean(6,:)]                   =   SensitivityAnalysis(Vars,Tcu_t,filename,'Sunlit Top Leaf Temperature',1);
[Vname,Vmean(7,:)]                   =   SensitivityAnalysis(Vars,Tcu_m,filename,'Sunlit Middle Leaf Temperature',1);
[Vname,Vmean(8,:)]                   =   SensitivityAnalysis(Vars,Tcu_b,filename,'Sunlit Bottom Leaf Temperature',1);

% average sensitivity
fid = fopen([outputdir,filename,'-Sensitivity.txt'],'w');
fprintf(fid,'Mean Sensitivity of \n');
fprintf(fid,'\t\tTsh\t\t\tTsu\t\t\tTch_t\t\tTch_m\t\tTch_b\t\tTcu_t\t\tTcu_m\t\tTcu_v\n');

for j=1:length(Vname)
    fprintf(fid,'%s:  ',Vname{j});
    fprintf(fid,'\t%4.3e',Vmean(:,j));
    fprintf(fid,'\n');
end
fclose(fid)

%% Perform Sensitivity analysis for MODIS bands
wlT                             =   spectral.wlS(iwls);
MODIS.mean = mean([MODIS.min;MODIS.max])

for j=1:length(MODIS.mean)
    fprintf(1,'%f - %f \n', MODIS.min(j), MODIS.max(j))
    MODISbandstr{j}             =   sprintf('MODIS @ %5.1fnm',MODIS.mean(j));
    
    wli                         =   MODIS.mean(j); linspace(MODIS.min(j),MODIS.max(j),5);    
    Lot_i                       =   interp1(wlT,double(Lot_)/100,wli);
    
       
    Lot_i                       =   reshape(Lot_i, Nsample); 
    [Vname,Vmean_m(j,:)]        =   SensitivityAnalysis(Vars,Lot_i,filename,MODISbandstr{j},1)

end