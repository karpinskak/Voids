clc
clear
clf

promienie=(1:1:87.75)*10^(-6); % zakres rozkladu prawdopodobienstwa dla 27 sierpnia
promienie_stat=[12.9-2*4.8,12.9-4.8,12.9,12.9+4.8,12.9+2*4.8]*10^(-6);
DIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Cloud_voids_theoretically/Stare_przydatne/';
load([DIR,'Data/Simulation_vortex_parameters.mat'],'vortex_param')

fDIR=[DIR '/Functions'];
addpath(DIR)
addpath(fDIR)
Const=Constants;
% punkty zawieraja wartosci bezwymiarowe
% punkty=[nr_promienia,promien, polozenie x, polozenie y, stabilnosc];

parpool('local',10)

parfor j=1:numel(vortex_param.A)
   clf
    punkty= [];
    for k=1:numel(promienie)
        [par_eq]=wylicz_param(Const,2,promienie(k),vortex_param.delta(j),vortex_param.teta(j),vortex_param.A(j),1,1);
        [points, stability] = eq_points(vortex_param.A(j),par_eq.Sv,par_eq.St);
        punkty_temp=[zeros(numel(points(:,1)),1)+k,(points(:,1).^2+points(:,2).^2).^(0.5),points,stability'];
        punkty=[punkty;punkty_temp];
    end
    %punkty_x=punkty(:,3)*vortex_param.delta(j);
    %punkty_y=punkty(:,4)*vortex_param.delta(j);
    wartosci_Sv=wykres(punkty(:,2),vortex_param.A(j));
     
    punkty_stat= [];
    for k=1:numel(promienie_stat)
        [par_eq]=wylicz_param(Const,2,promienie_stat(k),vortex_param.delta(j),vortex_param.teta(j),vortex_param.A(j),1,1);
        [points, stability] = eq_points(vortex_param.A(j),par_eq.Sv,par_eq.St);
        punkty_temp=[zeros(numel(points(:,1)),1)+k,(points(:,1).^2+points(:,2).^2).^(0.5),points,stability'];
        punkty_stat=[punkty_stat;punkty_temp];
    end
    
    wartosci_stat_Sv=wykres(punkty_stat(:,2),vortex_param.A(j));
    
    sztuczne_r=0.000001:0.1:1.1*max(max(punkty_stat(:,2)),15);
    cala_krzywa=wykres(sztuczne_r,vortex_param.A(j));
    
    punkty_fil=nonzeros(punkty(:,2).*((punkty(:,5)==1)));
    punkty_nfil=nonzeros(punkty(:,2).*((punkty(:,5)==2)));
    wartosci_fil=wykres(punkty_fil,vortex_param.A(j));
    wartosci_nfil=wykres(punkty_nfil,vortex_param.A(j));
    
    punkty_stat_fil=nonzeros(punkty_stat(:,2).*((punkty_stat(:,5)==1)));
    punkty_stat_nfil=nonzeros(punkty_stat(:,2).*((punkty_stat(:,5)==2)));
    wartosci_stat_fil=wykres(punkty_stat_fil,vortex_param.A(j));
    wartosci_stat_nfil=wykres(punkty_stat_nfil,vortex_param.A(j));
    
    figure(1)
    plot(sztuczne_r,cala_krzywa)
    hold on
    scatter(punkty_fil,wartosci_fil,20,'g','filled')
    hold on
    scatter(punkty_nfil,wartosci_nfil,15,'g')
    hold on
    scatter(punkty_stat_fil,wartosci_stat_fil,50,'r','filled')
    hold on
    scatter(punkty_stat_nfil,wartosci_stat_nfil,50,'r')
    hold off
    saveas(gcf,['/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Cloud_voids_theoretically/Stare_przydatne/Data/sim_eq_curves_plots/eq_curve_',num2str(j)],'fig')
    saveas(gcf,['/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Cloud_voids_theoretically/Stare_przydatne/Data/sim_eq_curves_plots/eq_curve_',num2str(j)],'png')

    %punkty_stat_x=punkty_stat(:,3)*vortex_param.delta(j);
    %punkty_stat_y=punkty_stat(:,4)*vortex_param.delta(j);
   eq_points_all{j}.punkty=punkty;
   eq_points_all{j}.punkty_stat=punkty_stat;
end

 save('Simulation_eq_points.mat', 'eq_points_all','promienie', 'promienie_stat','vortex_param')
 delete(gcp('nocreate'))
function wartosc=wykres(r,A)
wartosc=A*r.*sqrt(1+((1-exp(-r.^2/2))./(2*pi*A*r.^2)).^2);
end




                    
                    
                    