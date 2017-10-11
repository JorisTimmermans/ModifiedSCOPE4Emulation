%%
Tcu_mean                                    =   mean(mean(thermal.Tcu,1),2);

subplot(2,1,1,'nextplot','add')
plot(thermal.Tch,canopy.x,'Color',[0 0.5 0]) ;
plot(Tcu_mean(:),canopy.x,'Color',[0 1 0]) ;

plot(thermal.Ts(1),-1,'o','Color',[0.5 0 0]) ;
plot(thermal.Ts(2),-1,'o','Color',[1 0 0]) ;

plot(meteo.Ta,0,'o','Color',[0 0 1]) ;

legend('Shaded leaves','Sunlit leaves','Shaded soil', 'Sunlit soil','Air')  
subplot(2,1,2,'nextplot','add')
plot(spectral.wlS, rad.Lot_) ;
title(sprintf('%04.3f W/m2',rad.Lot))