% photosynthesis (A), fluorescence factor (F), and stomatal resistance (rcw), for shaded (1) and sunlit (h) leaves
% options.Fluorescence_model = 1;
%% Standard parameters
biochem_in.Type             = leafbio.Type;
if options.Fluorescence_model~=3
    biochem_in.Fluorescence_model = options.Fluorescence_model;
    biochem_in.p            = p;
    biochem_in.m            = leafbio.m;
    biochem_in.O            = meteo.Oa;
    biochem_in.Rdparam      = leafbio.Rdparam;
end

if options.Fluorescence_model==2    % specific for the v.Caemmerer-Magnani model
    b                       = @biochemical_MD12;
    biochem_in.Tyear        = leafbio.Tyear;
    biochem_in.beta         = leafbio.beta;
    biochem_in.qLs          = leafbio.qLs;
    biochem_in.NPQs         = leafbio.kNPQs;
    biochem_in.stressfactor = leafbio.stressfactor;
elseif options.Fluorescence_model==3
    b                       = @biochemical_JT16; % specific for Jarvis model (regressed to Berry - v.d. Tol model)
%     biochem_in.Tparams      = leafbio.Tparam;
%         biochem_in.stressfactor = SMCsf;    
else
    b                       = @biochemical;             % specific for Berry-v.d.Tol model
    biochem_in.tempcor      = options.apply_T_corr;
    biochem_in.Tparams      = leafbio.Tparam;
    biochem_in.stressfactor = SMCsf;    
end

if counter == 0
    biochem_in.A            = -999;                     % Net assimilation rate of the leaves
else 
    biochem_in.A            = Au;
end

%% for shaded leaves
% specify shaded leaves inputs
biochem_in.T                = Tch;                      % leaf temperature (required for all)
biochem_in.eb               = ech;                      % water vapor content (required for all)
biochem_in.Vcmo             = fVh.*leafbio.Vcmo;        % carboxylation (not required by JT)
biochem_in.Cs               = Cch;                      % CO2 concentration (not required for JT)
biochem_in.Q                = rad.Pnh*1E6;              % PAR(net radiation), Quantum Flux (required for all)

% Run Biochemical Model for shaded leaves
biochem_out                 = b(biochem_in);

% Retrieve outputs of the Model specific for shaded leaves
Ah                          = biochem_out.A;
Cih                         = biochem_out.Ci;
Fh                          = biochem_out.eta;
rcwh                        = biochem_out.rcw;
qEh                         = biochem_out.qE;       % vCaemmerer- Magnani does not generate this parameter (dummy value)

%% for sunlit leaves
biochem_in.T                = Tcu;
biochem_in.eb               = ecu;
biochem_in.Vcmo             = fVu.*leafbio.Vcmo;
biochem_in.Cs               = Ccu;
biochem_in.Q                = rad.Pnu*1E6;

% Run Biochemical Model for sunlit leaves
biochem_out                 = b(biochem_in);

% Retrieve outputs of the Model specific for sunlit leaves
Au                          = biochem_out.A;
Ciu                         = biochem_out.Ci;
Fu                          = biochem_out.eta;
rcwu                        = biochem_out.rcw;
qEu                         = biochem_out.qE;
    
% varnames    =   fieldnames(biochem_out);
% for i=1:length(varnames)
%     a = biochem_out.(varnames{i});
%     varnames{i}
%     if any(isnan(a(:)))
%         keyboard
%     end
% end
%% biophysical calculated rcw by simple fit 
% this was used in the 1st investigation (to check whether variables can be omitted)
% Pfitu                     =   [3.9132e-07   -6.8137e-05    0.0044];
% Pfith                     =   [4.9657e-07   -9.1193e-05    0.0047];
% 
% 
% mean_gcwh_fit             =   polyval(Pfitu,1:60);
% mean_gcwu_fit             =   polyval(Pfith,1:60);
% 
% gcwh_fit                  =   mean_gcwh_fit;
% gcwu_fit                  =   permute(meshgrid(mean_gcwu_fit,ones(13,1),ones(36,1)),[1 3 2]);
% 
% rcwh                      = 1./gcwh_fit';
% rcwu                      = 1./gcwu_fit;


    