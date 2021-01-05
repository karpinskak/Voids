DIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Cloud_voids_theoretically/Stare_przydatne/Results';
 %% find all files that are named "trajectories_all.fig"
 
file_list=dir(fullfile(DIR,'**/trajectories_all.fig'));

for j=1:numel(file_list)
    nazwa=[file_list(j).folder,'/trajectories_all'];
    openfig([file_list(j).folder,'/trajectories_all'])
    saveas(gcf,nazwa,'png');
    close all
end