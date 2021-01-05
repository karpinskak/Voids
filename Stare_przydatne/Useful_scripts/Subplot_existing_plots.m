clear
clc
close all

DIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Cloud_voids/';
spec1='Results/Pierwsze_kilka/A000036/Delta00078/R12_9/Th025/k6_41_l2_n10_tk6_poli_exper/';
spec2='Results/A00003/Delta001/R12_9/Th05/k10_00_l2_n10_tk8_poli_exper/';

nazwa_nR='n_r_Rvis_mean_plot';
load([DIR spec1 'Trajectories.mat'],'tglobal')
nazwa{1}=char(strcat(DIR, spec1, num2str(numel(tglobal)),'.fig'));
nazwa{5}=char(strcat(DIR, spec1, nazwa_nR,'.fig'));
%nazwa{3}=char(strcat(DIR, spec1, 'Fig', num2str(numel(tglobal)),'_3D.fig'));

load([DIR spec2 'Trajectories.mat'],'tglobal')
nazwa{2}=char(strcat(DIR, spec2, num2str(numel(tglobal)),'.fig'));
nazwa{6}=char(strcat(DIR, spec2, nazwa_nR,'.fig'));
%nazwa{3}=char(strcat(DIR, spec2, 'Fig', num2str(numel(tglobal)),'_3D.fig'));

% Load saved figures
for j=[1,2,5,6]%1:numel(nazwa)
figures{j}=hgload(nazwa{j});
end

% Prepare subplots
figure
for j=[1,2,5,6]%1:numel(nazwa)
h=subplot(3,2,j);
% Paste figures on the subplots
copyobj(allchild(get(figures{j},'CurrentAxes')),h);
end

% % Add legends
% l(1)=legend(h(1),'LegendForFirstFigure')
% l(2)=legend(h(2),'LegendForSecondFigure')