close all

SCOPE

iwl  = MODIS.nr>=31000 & MODIS.nr<=32000;


Tsh1 = thermal.Ts(2);
Tsu1 = thermal.Ts(1);
Tch1 = thermal.Tch;
Tcu1 = mean(mean(thermal.Tcu,1),2);
minV = min([min(Tsh1), min(Tsu1), min(Tch1), min(Tcu1), min(meteo.Ta)])*0.95;
maxV = max([max(Tsh1), max(Tsu1), max(Tch1), max(Tcu1), max(meteo.Ta)])*1.05;


figure('Position',[20 20 1024 800])
subplot(1,1,1,'Nextplot','add')
plot(Tsh1(:), 0 ,'o','color',[1.0,0,0])
plot(Tsu1(:), 0 ,'>','color',[0.5,0,0])

plot(Tch1(:),[1:60],'g','color',[0,0.5,0])
plot(Tcu1(:),[1:60],'g','color',[0,1.0,0])
plot(meteo.Ta,62 ,'db')
xlim([minV, maxV])
% subplot(2,1,2)
plot(MODIS.TEB(iwl),[70, 70],'mh')
xlim([minV, maxV])



Vcmo =  5;

SCOPE_translate_bck            
SCOPE_run

Tsh2 = thermal.Ts(2);
Tsu2 = thermal.Ts(1);
Tch2 = thermal.Tch;
Tcu2 = mean(mean(thermal.Tcu,1),2);

figure('Position',[20 20 1024 800])
subplot(1,1,1,'Nextplot','add')
plot(Tsh2(:), 0 ,'o','color',[1.0,0,0])
plot(Tsu2(:), 0 ,'>','color',[0.5,0,0])

plot(Tch2(:),[1:60],'g','color',[0,0.5,0])
plot(Tcu2(:),[1:60],'g','color',[0,1.0,0])
plot(meteo.Ta,62 ,'db')
xlim([minV, maxV])

% subplot(2,1,2)
plot(MODIS.TEB(iwl),[70, 70],'mh')
