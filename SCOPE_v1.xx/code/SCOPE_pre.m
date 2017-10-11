%% SCOPE.m (script)

%     SCOPE is a coupled radiative transfer and energy balance model
%     Copyright (C) 2015  Christiaan van der Tol
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%
clc
%clear all

%% 0. globals
global constants

%% 1. define constants
[constants] = define_constants();

%% 2. simulation options
path_of_code                = cd;
parameter_file              = {'input_data.xlsx'};  % USE XLSX file for cross-platform compatibility

try
    options                 = readStructFromExcel(['../' char(parameter_file)], 'options', 3, 1);
catch
    load('input_data.mat','options')
end

if options.simulation>2 || options.simulation<0, fprintf('\n simulation option should be between 0 and 2 \r'); return, end

%% 3. file names
try
    [dummy,X]                       = xlsread(['../' char(parameter_file)],'filenames');

    j = find(~strcmp(X(:,2),{''}));
    X = X(j,(1:end));
    F = struct('FileID',{'Simulation_Name','soil_file','leaf_file','atmos_file'...
        'Dataset_dir','t_file','year_file','Rin_file','Rli_file'...
        ,'p_file','Ta_file','ea_file','u_file','CO2_file','z_file','tts_file'...
        ,'LAI_file','hc_file','SMC_file','Vcmax_file','Cab_file'});
    for i = 1:length(F)
        k = find(strcmp(F(i).FileID,strtok(X(:,1))));
        if ~isempty(k)
            F(i).FileName = strtok(X(k,2));
            %if i==4, F(i).FileName = strtok(X(k,2:end)); end
        end
    end
catch
    
    load('input_data.mat','F')
end
%% 4. input data
try
    [N,X]                       = xlsread(['../' char(parameter_file)],'inputdata', '');    
    X                           = X(9:end,1);
    V                           = assignvarnames;
catch

    load('input_data.mat','V','X','N')
end    
    
options.Cca_function_of_Cab = 0;

for i = 1:length(V)
    j = find(strcmp(strtok(X(:,1)),V(i).Name));
    if isempty(j) || sum(~isnan(N(j,:)))<1
        if i==2
            fprintf(1,'%s %s %s \n','warning: input "', V(i).Name, '" not provided in input spreadsheet...');
            fprintf(1,'%s %s %s\n', 'I will use 0.25*Cab instead');
            options.Cca_function_of_Cab = 1;
        else

            if ~(options.simulation==1) && (i==30 || i==32)
                fprintf(1,'%s %s %s \n','warning: input "', V(i).Name, '" not provided in input spreadsheet...');
                fprintf(1,'%s %s %s\n', 'I will use the MODTRAN spectrum "',char(F(4).FileName(1)),'" as it is');
            else
                if (options.simulation == 1 || (options.simulation~=1 && (i<46 || i>50)))
                    fprintf(1,'%s %s %s \n','warning: input "', V(i).Name, '" not provided in input spreadsheet');
                    if (options.simulation ==1 && (i==1 ||i==9||i==22||i==23||i==54 || (i>29 && i<37)))
                        fprintf(1,'%s %s %s\n', 'I will look for the values in Dataset Directory "',char(F(5).FileName),'"');
                    else
                        if (i== 24 || i==25)
                            fprintf(1,'%s %s %s\n', 'will estimate it from LAI, CR, CD1, Psicor, and CSSOIL');
                            options.calc_zo = 1;
                        else
                            if (i>38 && i<44)
                                fprintf(1,'%s %s %s\n', 'will use the provided zo and d');
                                options.calc_zo = 0;
                            else
                                if ~(options.simulation ==1 && (i==30 ||i==32))
                                    fprintf(1,'%s \n', 'this input is required: SCOPE ends');
                                    return
                                else
                                    fprintf(1,'%s %s %s\n', '... no problem, I will find it in Dataset Directory "',char(F(5).FileName), '"');
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    if sum(~isnan(N(j,:)))<1
        V(i).Val            = -999;
    else
        V(i).Val            = N(j,~isnan(N(j,:)));
    end
end


%% 5. Declare paths
path_input      = '../../data/input/';          % path of all inputs

%% 6. Numerical parameters (iteration stops etc)
iter.maxit           = 300;                          %                   maximum number of iterations
iter.maxEBer         = 1;                            %[W m-2]            maximum accepted error in energy bal.
iter.Wc              = 1;                         %                   Weight coefficient for iterative calculation of Tc

%% 7. Load spectral data for leaf and soil
opticoef    = load([path_input,'fluspect_parameters/',char(F(3).FileName)]);  % file with leaf spectral parameters
rsfile      = load([path_input,'soil_spectrum/',char(F(2).FileName)]);        % file with soil reflectance spectra

% Optical coefficient data used by fluspect
optipar.nr    = opticoef(:,2);
optipar.Kab   = opticoef(:,3);
optipar.Kca   = opticoef(:,4);
optipar.Ks    = opticoef(:,5);
optipar.Kw    = opticoef(:,6);
optipar.Kdm   = opticoef(:,7);
optipar.phiI  = opticoef(:,9);
optipar.phiII = opticoef(:,10);
%optipar.GSV1  = opticoef(:,11); soil spectra, not used yet
%optipar.GSV2  = opticoef(:,12);
%optipar.GSV3  = opticoef(:,13);

%% 8. Load directional data from a file
if options.calc_directional
    anglesfile          = load([path_input,'directional/brdf_angles2.dat']); %     Multiple observation angles in case of BRDF calculation
    directional.tto     = anglesfile(:,1);              % [deg]             Observation zenith Angles for calcbrdf
    directional.psi     = anglesfile(:,2);              % [deg]             Observation zenith Angles for calcbrdf
    directional.noa     = length(directional.tto);      %                   Number of Observation Angles
end

%% 9. Define canopy structure
canopy.nlayers  = 60;
nl              = canopy.nlayers;
canopy.x        = (-1/nl : -1/nl : -1)';         % a column vector
canopy.xl       = [0; canopy.x];                 % add top level
canopy.nlincl   = 13;
canopy.nlazi    = 36;
canopy.litab    = [ 5:10:75 81:2:89 ]';   % a column, never change the angles unless 'ladgen' is also adapted
canopy.lazitab  = ( 5:10:355 );           % a row

%% 10. Define spectral regions
[spectral] = define_bands;

wlS  = spectral.wlS;    % SCOPE 1.40 definition
wlP  = spectral.wlP;    % PROSPECT (fluspect) range
wlT  = spectral.wlT;    % Thermal range
wlF  = spectral.wlF;    % Fluorescence range

I01  = find(wlS<min(wlF));   % zero-fill ranges for fluorescence
I02  = find(wlS>max(wlF));
N01  = length(I01);
N02  = length(I02);

nwlP = length(wlP);
nwlT = length(wlT);

nwlS = length(wlS);

spectral.IwlP = 1 : nwlP;
spectral.IwlT = nwlP+1 : nwlP+nwlT;
spectral.IwlF = (640:850)-399;

[rho,tau,rs] = deal(zeros(nwlP + nwlT,1));

% keyboardcl
%% 11. load time series data
if options.simulation == 1
    vi = ones(length(V),1);
    [soil,leafbio,canopy,meteo,angles,xyt]  = select_input(V,vi,canopy,options);
    [V,xyt,canopy]  = load_timeseries(V,leafbio,soil,canopy,meteo,constants,F,xyt,path_input,options);
else
    soil = struct;
end

%% 12. preparations
if options.simulation==1
    diff_tmin           =   abs(xyt.t-xyt.startDOY);
    diff_tmax           =   abs(xyt.t-xyt.endDOY);
    I_tmin              =   find(min(diff_tmin)==diff_tmin);
    I_tmax              =   find(min(diff_tmax)==diff_tmax);
    if options.soil_heat_method<2
        if (isempty(meteo.Ta) || meteo.Ta<-273), meteo.Ta = 20; end
        soil.Tsold = meteo.Ta*ones(12,2);
    end
end

nvars = length(V);
vmax = ones(nvars,1);
for i = 1:nvars
    vmax(i) = length(V(i).Val);
end
vmax([14,27],1) = 1; % these are Tparam and LADFb
vi      = ones(nvars,1);
switch options.simulation
    case 0, telmax = max(vmax);  [xyt.t,xyt.year]= deal(zeros(telmax,1));
    case 1, telmax = size(xyt.t,1);
    case 2, telmax  = prod(double(vmax)); [xyt.t,xyt.year]= deal(zeros(telmax,1));
end
[rad,thermal,fluxes] = initialize_output_structures(spectral);
atmfile     = [path_input 'radiationdata/' char(F(4).FileName(1))];
atmo.M      = aggreg(atmfile,spectral.SCOPEspec);

%% 13. create output files
create_output_files


%% preliminary stuff
Tb                          =   (-20:0.01:50);
dwl                         =   diff(spectral.wlS);                 %[nm]
DWL                         =   meshgrid(dwl,Tb);                   %[nm]

% MODIS.nr                    =   [  20    21    22   23    24    25   26    27    28  29   30    31    32     33     34     35     36]*1e3;
% MODIS.min                   =   [3.66 3.929 3.929 4.02 4.433 4.482 1.36 6.535 7.175 8.4 9.58 10.78 11.77 13.185 13.485 13.785 14.085]*1e3;
% MODIS.max                   =   [3.84 3.989 3.989 4.08 4.498 4.549 1.39 6.895 7.475 8.7 9.88 11.28 12.27 13.485 13.785 14.085 14.385]*1e3;

MODIS.nr                    =   [  20    21    22   23    24    25    27    28  29   30    31    32     33     34     35     36]*1e3;
MODIS.min                   =   [3.66 3.929 3.929 4.02 4.433 4.482 6.535 7.175 8.4 9.58 10.78 11.77 13.185 13.485 13.785 14.085]*1e3;
MODIS.max                   =   [3.84 3.989 3.989 4.08 4.498 4.549 6.895 7.475 8.7 9.88 11.28 12.27 13.485 13.785 14.085 14.385]*1e3;



bandwidth                   =   MODIS.max- MODIS.min;
ierror                      =   bandwidth<100;
MODIS.nr                    =   MODIS.nr(~ierror);
MODIS.min                   =   MODIS.min(~ierror);
MODIS.max                   =   MODIS.max(~ierror);

[~,iss]                      =   sort(MODIS.min);
MODIS.nr                    =   MODIS.nr(iss);
MODIS.min                   =   MODIS.min(iss);
MODIS.max                   =   MODIS.max(iss);

MODIS.LUT                   =   zeros(length(Tb), length(MODIS.nr)+1)*NaN;
MODIS.LUT(:,1)              =   Tb;
for iband = 1:length(MODIS.nr)
    iwl                     =   find(spectral.wlS>=MODIS.min(iband) & spectral.wlS<=MODIS.max(iband));
    wlT                     =   spectral.wlS(iwl);
    
    [WLT,TB]                =   meshgrid(wlT,Tb);
    
    
    Lb_hypers               =   Planck(WLT,TB + 273.15,1);           %[W m-2 str-1 um-1]
    dwlT                    =   DWL(:,iwl)*1e-3;            %[um]
    
    
    MODIS.dwlT{iband}       =   dwlT(1,:);
    
    MODIS_L                 =   sum(Lb_hypers.*dwlT,2);
    MODIS.LUT(:,iband+1)    =   MODIS_L;
    
    
    
%     hold on
%     plot(Tb-273.15,MODIS_L)
%     legend('1','2')
end

%%
% A=log(1./MODIS.LUT(:,2));
% B=MODIS.LUT(:,1);
% 
% [xData, yData] = prepareCurveData( A, B );
% 
% % Set up fittype and options.
% ft = fittype( 'smoothingspline' );
% 
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft );


%% define transmissitivities/reflectivities so that the output is consistent with MCD43B
% http://www.globalbedo.org/docs/GlobAlbedo_Albedo_ATBD_V4.12.pdf
% VIS (PAR) region [400 - 700nm] and NIR region [700 - 5000nm]
leaf_refPAR                                 =   0.0500;
leaf_transPAR                               =   0.0100;
soil_refPAR                                 =   0.1400;

leaf_refNIR                                 =   0.2529;
leaf_transNIR                               =   0.3316;
soil_refNIR                                 =   0.3487;

