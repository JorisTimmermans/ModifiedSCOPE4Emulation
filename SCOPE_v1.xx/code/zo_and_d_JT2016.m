function [z0m,d0] = zo_and_d_JT2016(canopy)
global constants
Cd                  =   canopy.Cd;
k                   =   constants.k;

LAI                 =   canopy.LAI;
hc                  =   canopy.hc;
dummy               =   zeros(size(LAI),'single');

%% U*/U(h)
c1                  =   0.320;                                                              % model constants (Massman 1997)
c2                  =   0.264;                                                              % model constants (Massman 1997)
c3                  =   15.1;                                                               % model constants (Massman 1997)
ust2u_h             =   c1 - c2 * exp(-c3 * Cd * LAI);                                      % Ratio of ustar and u(h) (Su. 2001 Eq. 8)
Cd_bs               =   2*ust2u_h.^2;                                                       % Bulk surface drag cofficient (Su 2001)

%% within-canopy wind speed profile extinction coefficient
n_ec                =   Cd * LAI ./ (Cd_bs);                                                % windspeed profile extinction coefficient (Su. 2002 Eq 7.)
In_ec               =   n_ec~=0;

d2h                 =   dummy;
d2h(In_ec)          =   1 - 1./(2*n_ec(In_ec)) .* (1 - exp(-2 * n_ec(In_ec)));              % Ratio of displacement height and canopy height (derived from Su 2002, eq 9)

%% Roughness length for momentum
z0m                 =   hc.*(1 - d2h) .* exp(-k * ust2u_h.^-1);                             % roughness heigth (eq 10, su 2001)

%% Displacement Height
d0                  =   d2h .* hc;                                                          % displacement height


%% Output
% canopy.zo           =   z0m;                                                                % roughness length for momentum of the canopy
% canopy.d            =   d0;                                                                 % displacement height of the canopy

