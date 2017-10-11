function [MODIS_TEB, LUT_error_min] =Simulate_MODIS_TEB(MODIS, rad, spectral)

[MODIS_TEB, LUT_error_min]  =   deal(zeros(length(MODIS.nr),1)*NaN);
for iband = 1:length(MODIS.nr)
    iwl                     =   find(spectral.wlS>=MODIS.min(iband) & spectral.wlS<=MODIS.max(iband));
%     wlT                     =   spectral.wlS(iwl);
    
    dwlT                    =   MODIS.dwlT{iband};
    
    Lband                   =   rad.Ltot_(iwl)';
%     Lband                   =   rad.Lot_(iwl)';
    
         
    MODIS_LL                 =   sum(Lband.*dwlT,2);        
    
    % find in LUT
    LUT_error               =   MODIS.LUT(:,iband+1)/pi-MODIS_LL;
    [LUT_error_min(iband),i]=   min(abs(LUT_error));
    MODIS_TEB(iband)        =   MODIS.LUT(i,1);
    
end