function [Vname,Vmean]=SensitivityAnalysis(Vars,Tsh,filename,Plotstring,close_option)


%%
outputdir                           =   '../output/Sensitivity_Analysis/';
if ~exist(outputdir,'dir')
    mkdir(outputdir)
end

%% Retrieve important variabels
names                               =   fieldnames(Vars);
Nvar                                =   length(names);

%% Define parameters for analysis
Q                                   =   [0.05, 0.25, 0.50, 0.75, 0.92];

%% Define indices
Nsample                             =   zeros(Nvar,1);
for iname = 1:Nvar;
    name                            =   names{iname};
    Nsample(iname)                  =   length(Vars.(name));
    I_bak.(name)                    =   1:Nsample(iname);
end

%% Loop over all variables
Nplots                              =   sum(Nsample>1);
Nc                                  =   3;
Nr                                  =   ceil(Nplots/Nc);
counter                             =   0;
S_q                                 =   cell(15,1);    

h1                                  =   figure('Position',[20 20 1024 800],'Visible','off');
h2                                  =   figure('Position',[20 20 1024 800],'Visible','off');

for iname = 1:Nvar
    name                            =   names{iname};
    x                               =   Vars.(name);
    n                               =   length(Vars.(name));
    if n>1  % only select variables that have been varied (otherwise no Sensitivity can be determined)
        counter                     =   counter +1;
        %% reshape data into Nx
        clear VV
        VV                          =   zeros(numel(Tsh)/n,n,'single');
        I                           =   I_bak;
        for i = 1:n
            I.(name)                =   i;

            V                       =   Tsh(I.Cab,I.Cw,I.Cca,I.Cdm,I.Cs,I.N,I.rho_thermal, I.tau_thermal, I.LAI, I.Ta, I.Rin, I.Rli, I.p, I.ea, I.u);
            VV(:,i)                 =   V(:);

        end
        VV(VV==0)                   =   NaN;

        %% Calculate Sensitivity as Jacobian
        dV                          =   diff(VV,[],2);    
        dx                          =   diff(x,[],2);
        xmean                       =   mean(cat(3,x(2:end),x(1:end-1)),3);
        dx                          =   repmat(dx,length(dV),1)*0+1;
        S                           =   dV./dx;
        
%         keyboard
        
        % quantiles
        V_q{iname}                  =   quantile(VV,Q,1);
        S_q{iname}                  =   quantile(abs(S),Q,1) + 1e-4;
%             keyboard
%         V_unc_lower_q90             =   S_q{iname}(1,:) - S_q{iname}(3,:);
%         V_unc_lower_q50             =   S_q{iname}(2,:) - S_q{iname}(3,:);
%         Vmean                       =   S_q{iname}(3,:) + 1e-3;
%         V_unc_higher_q50            =   S_q{iname}(4,:) - S_q{iname}(3,:);
%         V_unc_higher_q90            =   S_q{iname}(5,:) - S_q{iname}(3,:);

%         errorbar(xmean,[Vmean,V_unc_lower_q90,V_unc_higher_q90],'b')
%         errorbar(xmean,[Vmean,V_unc_lower_q50,V_unc_higher_q50],'r')
%         legend('q05-q95','q25-q75','mean','location','best')

        Vm  =   V_q{iname}(3,:);
%         Vl  =   V_q{iname}(3,:)-V_q{iname}(1,:);
        Vl2 =   V_q{iname}(3,:)-V_q{iname}(2,:);
%         Vu  =   V_q{iname}(5,:)-V_q{iname}(3,:);
        Vu2 =   V_q{iname}(4,:)-V_q{iname}(3,:);
        
        Vml  =   nanmin(VV)-Vm;
        Vmu  =   Vm-nanmax(VV);

        Sm  =   S_q{iname}(3,:);
%         Sl  =   S_q{iname}(3,:)-S_q{iname}(1,:);
        Sl2 =   S_q{iname}(3,:)-S_q{iname}(2,:);
%         Su  =   S_q{iname}(5,:)-S_q{iname}(3,:);
        Su2 =   S_q{iname}(4,:)-S_q{iname}(3,:);

        Sml  =   Sm - nanmin(abs(S));
        Smu  =   nanmax(abs(S))  - Sm;
        
        
        %% Pre Plotting stuff
        name2           =   name;
        name2(name2=='_')=   ' ';
                
        %% Plot Values
        h11     =   subplot(Nr,Nc,counter,'nextplot','add','Parent',h1);
        errorbar(x,Vm, Vml, Vmu,'Parent',h11,'b')
        errorbar(x,Vm, Vl2, Vu2,'Parent',h11,'r')
        plot(x,Vm,'g','Parent',h11,'linewidth',1)
        xlabel(h11,name2)
                
%         set(h2,'Visible','on')
        if counter==3
            legend(h11,'Extrema','Q25%-Q75%','median')
%             keyboard
        end
        if counter ==2
            title(h11,['Variation in ', Plotstring])
        end
        if mod(counter,Nc)==1
            ylabel(h11,['T [',char(176),'C]'])
        end
        
        
        
        %% Plot Sensitivity
        h(counter) = subplot(Nr,Nc,counter,'nextplot','add','yscale','log','Parent',h2);
        errorbar(xmean,Sm, Sml-1e-5, Smu,'Parent',h(counter))
        errorbar(xmean,Sm, Sl2-1e-5, Su2,'Parent',h(counter),'color','r')
        plot(xmean,Sm,'g','Parent',h(counter),'linewidth',1)
        
        xlabel(h(counter),name2)        
        if counter ==2
            title(h(counter),['Sensitivity of ', Plotstring])
        end
        if counter==3
            legend(h(counter),'Extrema','Q25%-Q75%','median')
%             keyboard
        end
        if mod(counter,Nc)==1
            ylabel(h(counter),['|\DeltaT| [' ,char(176),'C]'])
        end
        Vname{counter}  =   name;
        Vmean(counter)  =   nanmean(Sm);
        Vmax(counter)   =   nanmax(abs(S(:)));
        Vmin(counter)   =   nanmin(abs(S(:)));
        
        
    end
    
%     ylim([min(Vmin)*0.9 max(Vmax)*1.1])
    
end
linkaxes(h,'y')
ylim([1e-3, max(Vmax)*1.01])
set(h,'ytick',10.^[-3:1:5])
set(h,'ygrid','on','yminorgrid','off')


%%
print(h1,[outputdir,filename,'-Variation in ',Plotstring,'.png'],'-dpng','-r300')
print(h2,[outputdir,filename,'-Sensitivity in ',Plotstring,'.png'],'-dpng','-r300')

if close_option==1
    close all;
else
    set(h1,'Visible','on')
    set(h2,'Visible','on')
end
%% Output data
fid = fopen([outputdir,filename,'-Sensitivity in ',Plotstring,'.txt'],'w');
fprintf(fid,'Mean Sensitivity in %s\n',Plotstring);
for j=1:length(Vname)
    fprintf(fid,'%s: %4.3e\n',Vname{j},Vmean(j));
end
fclose(fid);
% edit([outputdir,'Sensitivity in ',Plotstring,'.txt'])

pause(1)


