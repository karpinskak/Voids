% Load constants and functions
addpath(DIR)
Const=Constants;
fDIR=[DIR 'Functions'];
addpath(fDIR)

% plot parameters
skala=100;
fsize=16;
width=16;
height=30;

figure(1)
set(gcf,'Position', [640, 300, 2*560, 2*420 ],...
    'paperunits','centimeters',...
    'Color',[1 1 1],...
    'papersize',[width,height],...
    'InvertHardCopy','off')
box on
[Delta,AA]=meshgrid(delta,A);
if size(R_2{1})~=size(Delta)
   for j=1:numel(teta)
       R_2{j}=R_2{j}';
   end
end

tiledlayout(3,numel(teta))

h(1)=nexttile;
hp(1)=pcolor(Delta*skala,AA,R_1{1});
ylabel('$A$','FontSize',fsize+2,'interpreter','latex')
grid on
h(2)=nexttile;
hp(4)=pcolor(Delta*skala,AA,R_1{2});
grid on
h(5)=nexttile;
hp(7)=pcolor(Delta*skala,AA,R_1{3});
ylabel('$A$','FontSize',fsize+2,'interpreter','latex')
grid on
h(3)=nexttile;
hp(2)=pcolor(Delta*skala,AA,R_2{1});
ylabel('$A$','FontSize',fsize+2,'interpreter','latex')
grid on
h(4)=nexttile;
hp(5)=pcolor(Delta*skala,AA,R_2{2});
grid on
h(6)=nexttile;
hp(8)=pcolor(Delta*skala,AA,R_2{3});
grid on
 d(1)=nexttile;
hp(3)=pcolor(Delta*skala,AA,R_2{1}-R_1{1});
xlabel('$\delta$~$[cm]$','FontSize',fsize+2,'interpreter','latex')
ylabel('$A$','FontSize',fsize+2,'interpreter','latex')
grid on
d(2)=nexttile;
hp(6)=pcolor(Delta*skala,AA,R_2{2}-R_1{2});
xlabel('$\delta$~$[cm]$','FontSize',fsize+2,'interpreter','latex')
grid on
d(3)=nexttile;
hp(9)=pcolor(Delta*skala,AA,R_2{3}-R_1{3});
xlabel('$\delta$~$[cm]$','FontSize',fsize+2,'interpreter','latex')
grid on

annotation('textbox',[1 0.5 0.3 0.2],'String','$\theta=0$')
%annotation('textbox',[],'$\theta=\frac{pi}{4}$','FontSize',fsize+2,'interpreter','latex')
%annotation('textbox',[],'\theta=\frac{pi}{2}$','FontSize',fsize+2,'interpreter','latex')
set(h,'FontSize',fsize)
set(d,'FontSize',fsize)
box on
set(hp,'EdgeColor','none')
set(h,'colormap',jet,'clim',[0 50])
cbh=colorbar(h(numel(h)-1));
title(cbh,'$R_<$~$[\mu m]$','FontSize',fsize,'interpreter','latex')
set(h,'colormap',jet,'clim',[0 50])
cbh=colorbar(h(numel(h)));
title(cbh,'$R_>$~$[\mu m]$','FontSize',fsize,'interpreter','latex')
set(d,'colormap',othercolor('BuDRd_12'),'clim',[0 50])
cbd=colorbar(d(numel(d)));
title(cbd,'$\Delta R$~$[\mu m]$','FontSize',fsize,'interpreter','latex')
