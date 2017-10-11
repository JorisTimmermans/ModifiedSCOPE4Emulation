function [spectral] = define_bands

    % Define spectral regions for SCOPE v_1.40
    % All spectral regions are defined here as row vectors
    % WV Jan. 2013
    
    % 3 spectral regions for SCOPE
    res_VIS         =      1;
    res_TIR         =    100;
    res_FTIR        =   1000;
    reg1            =   400 :   res_VIS  :  2400;
    reg2            =  2500 :   res_TIR  : 15000;
    reg3            = 16000 :   res_FTIR : 50000;
    
    
    spectral.wlS  = [reg1 reg2 reg3];
        
    % Other spectral (sub)regions    
    spectral.wlP   = reg1;                            % PROSPECT data range
    spectral.wlE   = 400:1:750;                       % excitation in E-F matrix
    spectral.wlF   = 640:1:850;                       % chlorophyll fluorescence in E-F matrix
    spectral.wlO   = reg1;                            % optical part
    spectral.wlT   = [reg2 reg3];                     % thermal part
    
    wlS            = spectral.wlS;
    
    spectral.IwlPAR =   wlS>=400 & wlS<=700;
    spectral.wlPAR  =   wlS(spectral.IwlPAR);           % PAR range
    spectral.IwlNIR =   wlS>=700 & wlS<=5000;           % NIR region (as defined by http://www.globalbedo.org/docs/GlobAlbedo_Albedo_ATBD_V4.12.pdf)
    
    % Data used by aggreg routine to read in MODTRAN data
    
    spectral.SCOPEspec.nreg = 3;
    spectral.SCOPEspec.start = [ 400  2500  16000];
    spectral.SCOPEspec.end   = [2400 15000  50000];
    spectral.SCOPEspec.res   = [res_VIS   res_TIR   res_FTIR];
    
    %% subset bands (not working yet with the rest of SCOPE
%     spectral.wlS                =   subset_spectral(spectral.wlS);
%     spectral.wlT                =   subset_spectral(spectral.wlT);
%     spectral.wlT_sub            =   subset_spectral(spectral.wlT);
end


