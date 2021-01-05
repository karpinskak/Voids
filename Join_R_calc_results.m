clc
clear
close all

src_folder='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Cloud_voids_theoretically/Results_Void_R_range/';
nzw1='_1_teta';
nzw2='_3teta';
dst_folder='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Cloud_voids_theoretically/';

x_ozn={'1','2','s'};
addpath(src_folder)


%% R1
load([src_folder,'R1_void',nzw1,'.mat'])
teta_1=teta;
R_1_1=R_1;
load([src_folder,'R1_void',nzw2,'.mat'])
R_1_2=R_1;
teta_2=teta;
clear R_1 teta
R_1=[R_1_2,R_1_1];
teta=[teta_2,teta_1];
clear R_1_1 R_1_2 teta_2 teta_1
save([dst_folder,'R1_void','.mat'],'A','delta','teta','R_1')

%% R2
clear R_1 teta 
load([src_folder,'R2_void',nzw1,'.mat'])
teta_1=teta;
R_2_1=R_2;
load([src_folder,'R2_void',nzw2,'.mat'])
R_2_2=R_2;
teta_2=teta;
clear R_2 teta
R_2=[R_2_2,R_2_1];
teta=[teta_2,teta_1];
clear R_2_1 R_2_2 teta_2 teta_1
save([dst_folder,'R2_void','.mat'],'A','delta','teta','R_2')

%% Rs
clear R_2 teta 
load([src_folder,'Rs_void',nzw1,'.mat'])
teta_1=teta;
R_s_1=R_s;
load([src_folder,'Rs_void',nzw2,'.mat'])
R_s_2=R_s;
teta_2=teta;
clear R_s teta
R_s=[R_s_2,R_s_1];
teta=[teta_2,teta_1];
clear R_s_2 R_s_1 teta_2 teta_1
save([dst_folder,'Rs_void','.mat'],'A','delta','teta','R_s')
