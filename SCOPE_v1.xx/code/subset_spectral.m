function [wl_sub] = subset_spectral(wl)
%% define spectral characteristics of the Sensors
Sensor.MODIS_TERRA.band_min   =   [620 841 459 545 1230 1628 2105 405 438 483 526 546 662 673 743 862 890 931 915 ... 
                             3660 3929 3929 4020 4433 4482 1360 6535 7175 8400 9580 10780 11770 13185 13485 13785 14085];
Sensor.MODIS_TERRA.band_max   =   [670 876 479 565 1250 1652 2155 420 448 493 536 556 672 683 753 877 920 941 965 ...
                             3840 3989 3989 4080 4498 4549 1390 6895 7475 8700 9880 11280 12270 13485 13785 14085 14385];
Sensor.MODIS_TERRA.band_nr     =   1:length(Sensor.MODIS_TERRA.band_min);

% derived from S3-TN-ESA-PL-316_SLSTR_Wavelengths_and_Irradiances_2012, derived by specifying SPRmin=0.01 
Sensor.S3_SLSTR.band_min      =   [543   647   850  1361  1569  2214  3365 10105 11165]/1000; %derived from 
Sensor.S3_SLSTR.band_max      =   [568   673   880  1391  1654  2289  4105 11495 13105]/1000;
Sensor.S3_SLSTR.band_nr     =   1:length(Sensor.S3_SLSTR.band_min);

Sensor.S3_OLCI.band_min       =   [391   405   435   484   504   553   613   658   668   676   702   748   758   761   765   770   853   878   893   928   998];
Sensor.S3_OLCI.band_max       =   [409   419   449   497   517   567   627   672   679   687   715   759   764   768   770   788   877   892   907   952  1041];
Sensor.S3_OLCI.band_nr         =   1:length(Sensor.S3_OLCI.band_min);

%% determine if wavelength is in one of the bands
sensornames             =   fieldnames(Sensor);
for iwl=1:length(wl)
    for isensor=1:length(sensornames)
        sensorname      =   sensornames{isensor};
        In(isensor,iwl)     =   any(wl(iwl)>=Sensor.(sensorname).band_min & wl(iwl)<=Sensor.(sensorname).band_max);
        
    end
end

In                         =   any(In,1);
wl_sub                     =   wl(In);

end