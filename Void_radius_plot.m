clear
close all
clc
delete(gcp('nocreate'))

DIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Cloud_voids_theoretically/';

% Load constants and functions
addpath(DIR)
Const=Constants;
fDIR=[DIR 'Functions/'];
addpath(fDIR)
load([DIR,'Stare_przydatne/Data/Simulation_vortex_parameters.mat'],'vortex_param')

% plot parameters
skala=100;
fsize=16;
width=16;
height=30;
v1=[0.5,2,5]; % wektor wartosci konturow
linest={'-','--',':'};
kolory=[[0;0.5;1],[1;0;0],[0.1;0.8;0]];%,[0;0;1],[1;230/255;0]];
kolory2=[[0.5;0.5;1],[1;0.5;0.5],[0.5;1;0.5]];
sim_points=1;
point_size=30;

%% Input data and calculate param=St/A;

npool=13;
nazwa='orbit_radii_data';
delta=[0.1:0.005:1.005]*10^(-2);
A=0.0001:2*10^(-4):0.03001;
%delta=[0.1:0.02:1]*10^(-2);
%A=0.0001:(1*10^(-4)):0.002;
R=[3,13,23]*10^(-6);
[Delta,AA]=meshgrid(delta,A);
skala=100;


ST=zeros(size(AA));
if isfile([DIR nazwa '.mat'])==1
    load([DIR nazwa])
else
    tic
    parpool('local',npool)
    for j=1:numel(R)
        for k=1:size(AA,1)
            for l=1:size(AA,2)
                [par]=wylicz_param(Const,2,R(j),Delta(k,l),0,AA(k,l),0,0);
                ST(k,l)=par.St;
            end
        end
        STT{j}=ST;
        clear par
        param=ST./AA;
        r0=zeros(size(AA));
        for k=1:size(param,1)
            r0_temp=zeros(1,size(param,2))';
            parfor l=1:size(param,2)
                ptemp=param(k,l);
                if ptemp>Const.St_A_cr
                    r0_temp(l)=Row_orb_bur(1/ptemp).*Delta(k,l);
                else
                    r0_temp(l)=NaN;
                end
            end
            r0(k,:)=r0_temp;
        end
        orbit_radii{j}=r0;
    end
    toc
    delete(gcp('nocreate'))
    save([DIR, nazwa],'STT','AA','A','Delta','delta','R','orbit_radii')
end
tau_p=2*R.^2*Const.ro_p/(9*Const.ro_a*Const.nu);

%% plot
figure(1)
set(gcf,'Position', [640, 300, 2*560, 2*420 ],...
    'paperunits','centimeters',...
    'Color',[1 1 1],...
    'papersize',[width,height],...
    'InvertHardCopy','off')
box on

for p=1:numel(R)
    j=numel(R)-p+1;
    Acrit=(tau_p(j)*Const.nu./delta.^2).^(0.5)/(4*pi);
    area(delta*skala,Acrit,'FaceColor',kolory2(:,j))
    alpha 0.5
    hold on
    
    r0temp=orbit_radii{j};
    Rtemp(p)=R(j);
    Delta(r0temp==0)=NaN;
    AA(r0temp==0)=NaN; % rysowanie tylko tych przypadków, dla których istnieje orbita stacjonarna
    C=r0temp*skala;
    
    [Cx,h]=contour(Delta*skala,AA,C,v1,'Color',kolory(:,j),'ShowText','off','Linewidth',3,'LineStyle','-');
    clabel(Cx,h,'FontSize',13,'FontName','Helvetica')
    hold on
end

legend({'R=23um','','R=13um','','R=3um',''})%,'simulations'},'FontSize',19,'FontName','Helvetica')
hold on
if sim_points==1
    lab=1:numel(vortex_param.delta);
    labelki=cellstr(num2str(lab'));
    scatter(vortex_param.delta*skala,vortex_param.A,point_size,'k','filled')
    text(1.01*vortex_param.delta*skala,1.01*vortex_param.A,labelki',30)
end
%xlim([10^(-1) 10^(1)])
%ylim([10^(-6) 10^(-2)])
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'fontsize',fsize+2)
xlabel('$\delta$~ $[cm]$','interpreter','latex','FontSize',40)
ylabel('$A$','interpreter','latex','FontSize',40)


% x=[0.5,0.2,0.5];
% y=[0.00036,0.0048,0.0076]; % symulacja 1) dziura czesciowa 2) dziura wyrazna 3) brak dziury
% scatter(x,y,70,'filled','MarkerFaceColor','black')
% tekst={'1)','2)','3)'};
% for s=1:3
% text(0.96*x(s),0.8*y(s),tekst(s),'FontSize',16,'FontWeight','bold','FontName','Helvetica')
% hold on
% end
