clear
close all
clc
delete(gcp('nocreate'))

DIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Cloud_voids/';

% Load constants and functions
addpath(DIR)
Const=Constants;
fDIR=[DIR 'Functions/'];
addpath(fDIR)
nazwar1='R1_void';
nazwar2='R2_void_x';
nazwar3='Rs_void';

% plot parameters
skala=100;
fsize=16;
width=16;
height=30;

if isfile([DIR nazwar1, '.mat']) && isfile([DIR nazwar2, '.mat']) && isfile([DIR nazwar3, '.mat'])
    load([DIR nazwar1,'.mat'])
    load([DIR nazwar2,'.mat'])
    load([DIR nazwar3,'.mat'])
else
    error('Lacking data file.')
end

[Delta,AA]=meshgrid(delta,A);

R_l=cell(numel(teta),1);
R_r=cell(numel(teta),1);

for j=1:numel(teta)
    R_ltemp=R_1{j};
    R_ltemp(R_1{j}>=R_s{j})=NaN;
    R_l{j}=R_ltemp;
    clear R_ltemp
    R_rtemp=R_2{j}.*(R_2{j}<R_s{j})+R_s{j}.*(R_2{j}>=R_s{j});
    R_r{j}=R_rtemp;
end

figure(1)
set(gcf,'Position', [640, 300, 2*560, 2*420 ],...
    'paperunits','centimeters',...
    'Color',[1 1 1],...
    'papersize',[width,height],...
    'InvertHardCopy','off')
box on

tiledlayout(3,numel(teta))

for j=1:numel(teta)
h(j)=nexttile;
hp(j)=pcolor(Delta*skala,AA,R_l{j});
if mod(j,3)==1
ylabel('$A$','FontSize',fsize+2,'interpreter','latex')
end
grid on
title(['$\theta=\pi/$',num2str(floor(1/(teta(j)/pi)))],'FontSize',fsize+2,'interpreter','latex')
end

for j=1:numel(teta)
h(j+numel(teta))=nexttile;
hp(j+numel(teta))=pcolor(Delta*skala,AA,R_r{j});
if mod(j,numel(teta))==1
ylabel('$A$','FontSize',fsize+2,'interpreter','latex')
end
grid on
end

for j=1:numel(teta)
d(j)=nexttile;
interv=R_r{j}-R_l{j};
interv(interv<0)=NaN;
hp(j+2*numel(teta))=pcolor(Delta*skala,AA,interv);
if mod(j,numel(teta))==1
ylabel('$A$','FontSize',fsize+2,'interpreter','latex')
end
grid on
xlabel('$\delta$~$[cm]$','FontSize',fsize+2,'interpreter','latex')

end

set(h,'FontSize',fsize-2)
set(d,'FontSize',fsize-2)
set(hp,'EdgeColor','none')
set(h,'yscale','log')
set(h,'xscale','log')
set(d,'yscale','log')
set(d,'xscale','log')
box on

set(h,'colormap',jet,'clim',[0 50])
cbh=colorbar(h(numel(teta)));
cbh.Label.String='$R_<$~$[\mu m]$';
cbh.Label.Interpreter='latex';
cbh.Label.FontSize=fsize;

cbh=colorbar(h(2*numel(teta)));
cbh.Label.String='$R_>$~$[\mu m]$';
cbh.Label.Interpreter='latex';
cbh.Label.FontSize=fsize;


set(d,'colormap',jet,'clim',[0 50])
cbd=colorbar(d(numel(teta)));
cbd.Label.String='$\Delta R$~$[\mu m]$';
cbd.Label.Interpreter='latex';
cbd.Label.FontSize=fsize;
