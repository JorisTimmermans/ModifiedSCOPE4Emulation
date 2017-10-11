function [biochem_out] = biochemical_JT16(biochem_in)
% Date: 	11 Oct 2017
% updated: 	04 04 2017 (added Vcmo-function)
%   
% Authors: 	Joris Timmermans 
% Sources: 	Jarvis 1976, Stewart 1988, White 2000, Damour 2010
%
% This function calculates:
%    - stomatal resistance of a leaf or needle (s m-1), 
%
% Information:
% This function calculates the stomatal resistance on basis of the Jarvis model. The baseline functions
% within the model are based on function-fitting of the original biochemical.m. This original module was 
% run for more than 10.000 scenarios. Afterwards, the convex hull of the point cloud was determined for
% each variable. Afterwards the specific baseline functions was fitted to these convex hulls.
%
% Usage:
% biochem_out = biochemical(biochem_in)
% the function was tested for Matlab 7.2.0.232 (R2006a)
%
%
% Input (units are important):
% All inputs are optional. structure 'biochem_in' with the following elements:
% Q         % [umol photons m-2 s-1]net radiation, PAR
% T         % [K]             leaf temperature
% eb        % [hPa]                 intial estimate of the vapour pressure in leaf boundary layer
% p         % [hPa]                 air pressure
% Type      % []                    text parameter, either 'C3' for C3 or any other text for C4
% stressfactor []                   optional input: stress factor to reduce Vcmax (for
%                                   example soil moisture, leaf age). Default value = 1.
%
% Note: always use the prescribed units. Temperature can be either oC or K
% Note: input can be single numbers, vectors, or n-dimensional
% matrices
%
% Output:
% structure 'biochem_out' with the following elements:
% A         % [umol/m2/s]           net assimilation rate of the leaves
% Cs        % [umol/m3]             CO2 concentration in the boundary layer
% eta0      % []                    fluorescence as fraction of dark
%                                   ...adapted (fs/fo)
% rcw       % [s m-1]               stomatal resistance
% qE        % []                    non photochemical quenching
% fs        % []                    fluorescence as fraction of PAR
% Ci        % [umol/m3]             internal CO2 concentration
% Kn        % []                    rate constant for excess heat
% fo        % []                    dark adapted fluorescence (fraction of aPAR)
% fm        % []                    light saturated fluorescence (fraction of aPAR)
% qQ        % []                    photochemical quenching
% Vcmax     % [umol/m2/s]           carboxylation capacity after
%                                   ... temperature correction 

if nargin<1
    biochem_in.Type         =   'C3';
end
Type                        =   biochem_in.Type;

%% Quantum Flux, defined as net PAR
if isfield(biochem_in,'Q')
    % Coefficients (with 95% confidence bounds):
    switch Type
        case 'C3'
            
            a               =   0.9936; %  (0.7742, 1.213)
            b               =   -1.323e-05; %  (-0.0001396, 0.0001132)
            c               =   -0.981; %  (-1.191, -0.7715)
            d               =   -0.003492; %  (-0.004765, -0.002219)
        otherwise
            % Coefficients (with 95% confidence bounds):
            a               =   0.6732; %  (0.6732, 0.6732)
            b               =   0.0002111; %  (0.0002111, 0.0002111)
            c               =   -0.6885; %  (-0.6885, -0.6885)
            d               =   -0.00332; %  (-0.00332, -0.00332)            
    end
    thresh                  =   2400;
    
    % model
    Q                       =   biochem_in.Q;
    kQ                      =   a*exp(b*Q) + c*exp(d*Q);
    
    diff                    =   Q-thresh;
    [~,ix]                  =   min(abs(diff(:)));
    
    kQ(Q>thresh)            =   kQ(ix);
    kQ                      =   max(kQ,0);
    
else
    kQ                      =   1;
end

%% Vcmo
if isfield(biochem_in,'Vcmo')
    %Coefficients (with 95% confidence bounds):
    switch Type
        case 'C3'
            p1              =   -7.44e-05; % (-9.361e-05, -5.519e-05)
            p2              =     0.01905; % (0.01741, 0.0207)
            p3              =    -0.08012; % (-0.1034, -0.05687)
    
        otherwise
    %         Coefficients (with 95% confidence bounds):
            p1              =  -0.0001523; % (-0.0002288, -7.578e-05)
            p2              =     0.02647; % (0.01806, 0.03488)
            p3              =     -0.1479; % (-0.3056, 0.009855)

    end

    % model
    Vcmo                    =   biochem_in.Vcmo;
    kVcmo                   =   p1*Vcmo.^2 + p2*Vcmo + p3;
    kVcmo                   =   min(kVcmo,1);
    kVcmo                   =   max(kVcmo,0);

    
else
    kVcmo                   =   1;
end
%% Temperature
if isfield(biochem_in,'T')
    switch Type
        case 'C3'
            % Coefficients (with 95% confidence bounds):
            a1              =	0.9501;%  (0.8922, 1.008)
            b1              =   0.05126;%  (0.04885, 0.05366)
            c1              =   -0.9286;%  (-1.658, -0.1989)
            a2              =   0.1134;%  (0.06255, 0.1642)
            b2              =   0.1765;%  (0.146, 0.2071)
            c2              =   11.14;%  (1.923, 20.36)
            
            Tthresh         =   [261,322];
            maxkT           =   1.0527;
        otherwise
            % Coefficients (with 95% confidence bounds):
            a1              =   8.132; % (-1737, 1754)
            b1              =   0.09273; % (-0.3389, 0.5243)
            c1              =   -6.87; % (-127, 113.3)
            a2              =   7.63; % (-1738, 1753)
            b2              =   0.09678; % (-0.3729, 0.5664)
            c2              =   51.69; % (-80, 183.4)
            
            Tthresh         =   [274,323];
            maxkT           =   1.0178;
    end

    % model
    T                       =   biochem_in.T + 273.15;          %temperature in K
    kT                      =   a1*sin(b1*T+c1) + a2*sin(b2*T+c2);
    
    
    ithresh                 =   T<Tthresh(1) | T>Tthresh(2);
    kT(ithresh)             =   1e-4;
    kT                      =   kT/maxkT;
    
else
    kT                      =   1;
end

%% water vapor
if isfield(biochem_in,'eb')
    % Coefficients (with 95% confidence bounds):
    switch Type
        case 'C3'            
            p1              =   6.83e-06;%  (2.675e-06, 1.098e-05)
            p2              =   -0.001352;%  (-0.001996, -0.0007075)
            p3              =   0.06621;%  (0.04033, 0.09209)
            p4              =   0.09196;%  (-0.1282, 0.3121)% Linear model Poly3:
            
            maxkE           =   1.0502; 
            eb_thresh       =   99;
        otherwise
            p1              =   1.892e-06; % (8.83e-07, 2.902e-06)
            p2              =   -0.0006115; % (-0.0007867, -0.0004363)
            p3              =   0.04605; % (0.03807, 0.05402)
            p4              =   -0.04158; % (-0.139, 0.05582)
            
            maxkE           =   1.9430;
            eb_thresh       =   110;
    end
    
    % model    
    eb                      =   biochem_in.eb;
    kE                      =   p1*eb.^3 + p2*eb.^2 + p3*eb + p4;
    kE(eb>eb_thresh)        =   0;
    
    kE                      =   kE/maxkE;
else
    kE                      =   1;
end

%% Pressure
if isfield(biochem_in,'p')
    % Coefficients (with 95% confidence bounds):
    switch Type
        case 'C3'            
            p1              =   -0.0006403;%  (-0.002096, 0.0008158)
            p2              =   1.512;%  (0.1918, 2.832)
        otherwise
            p1              =   -0.001962; %  (-0.01343, 0.009506)
            p2              =   2.566; %   (-7.604, 12.74)
    end

    % model:
    P                       =   biochem_in.p;
    kP                      =   p1*P + p2;
    kP                      =   min(kP,1);
    kP                      =   max(kP,0.5);
else
    kP                      =   1;
end

%% stress factor
if isfield(biochem_in,'stressfactor')
    switch Type
        case 'C3'
            a               =   NaN;
            b               =   NaN;
        otherwise
            a               =   NaN;
            b               =   NaN;
    end
    % model:
    s                       =   biochem_in.stressfactor;
    ks                      =   s*a + b;
    keyboard
else
    ks                      =   1;
end

%% Final
gsmax                       =   0.0079;
gs                          =   gsmax .* kQ .* kT .* kE .* kP .* kVcmo .*ks;
rcw                         =   1./gs;

ierror                      =   rcw<0 | rcw>0.625*1E6 | isnan(rcw);
rcw(ierror)              	=   0.625*1E6;

biochem_out.rcw             =   rcw;

dummy                       =   rcw*-9999;
fieldnames                  =   {'A','Cs','eta','qE','fs','Ci','Kn','fo','fm','qQ','Vcmax'};
for ivar =1:length(fieldnames)
    varname                     =   fieldnames{ivar};
    biochem_out.(varname)       =   dummy;
end

% if std(rcw(:))>200
%     keyboard
% end
