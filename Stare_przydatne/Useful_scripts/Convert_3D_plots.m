clear
close all
clc
delete(gcp('nocreate'))
tic

%% Data
DIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Cloud_voids/';

spec_dir={...
    %'Results/A00014/Delta001/R12_9/Th05/k5_l2_n10_tk8_poli_exper/',...
    %'Results/A00003/Delta001/R12_9/Th075/k5_l2_n10_tk8_poli_exper/',...
    'Results/A0003/Delta0005/R12_9/Th075/k10_l2_n10_tk6_poli_exper/',
    };
skala=100;
D=0.05;
Z=0.1;
for j=1
    loadDIR=char(strcat(DIR,spec_dir(j)));
    load([loadDIR,'Trajectories.mat'],'tglobal')
    
    for k=2541:20:numel(tglobal)
        
        nazwa=char(strcat(loadDIR, 'Animation_3D/Fig', num2str(k), '_3D.fig'));
        h=openfig(nazwa,'invisible');
        
        xmin=-1.1*D*skala;
        xmax=1.1*D*skala;
        ymin=-1.1*D*skala;
        ymax=1.1*D*skala;
        zmin=-1.1*Z*skala;
        zmax=1.1*Z*skala;
        
        set(gca,'XTick',linspace(ceil(xmin),floor(xmax),3));
        set(gca,'YTick',linspace(ceil(ymin),floor(ymax),3));
        set(gca,'ZTick',linspace(ceil(zmin),floor(zmax),12));
        
        [az,el] = view;
        view(az,el*1.0001)
        nazwa1=char(strcat(loadDIR, 'Animation_3D/Fig', num2str(k), '_3D'));
        saveas(h,nazwa1,'fig')
        saveas(h,nazwa1,'png')
        close(h)
        
    end
end