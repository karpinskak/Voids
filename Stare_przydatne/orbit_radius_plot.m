clear
clc
close all

% Parameters
v1=[0.5,1.5,3]; % wektor wartosci konturow w cm

% Plot parameters
width=16;
height=20;
linest={'-','--',':'};
kolory=[[0.4;0.8;1],[1;0.7;0.7],[0.5;1;0.5]];
delta_min=0.001;
delta_max=0.0175;
Amin=5*10^(-6);
Amax=0.01;
skala=100;

%% Load data
Con=Constants;
DIR=pwd;
load([DIR '/Dane/' 'promienie_orbit_A_gamma_Dzien1.mat'])
Rr=promienie_orbit.R;
gamma=promienie_orbit.gamma;
G=promienie_orbit.G;
Aa=promienie_orbit.A;
AA=promienie_orbit.AA;
r0=promienie_orbit.r0;
tau_a=promienie_orbit.taua;

% obciecie rysowanego zakresu
gmin=2*Con.nu./delta_max^2;
gmax=2*Con.nu./delta_min^2;
[~,I_gmin]=min(abs(gamma-gmin)); %pozycja minimalnej wartosci gammy po obcieciu
[~,I_gmax]=min(abs(gamma-gmax));%;
I_Amin=find(Aa==Amin);%;
I_Amax=find(Aa==Amax);% %pozycja maksymalnej wartosci Lw po obcieciu
gamma=gamma(I_gmin:I_gmax);
delta=(2*Con.nu./gamma).^(1/2);
Aa=Aa(I_Amin:I_Amax);


figure(1)
set(gcf,'Position', [640, 300, 2*560, 2*420 ])

for p=1:numel(r0)
    j=numel(r0)-p+1; % 3,2,1
    Acrit=(tau_a(j)*gamma/2).^(0.5)/(4*pi);
    plocik=area(fliplr(delta*skala),fliplr(Acrit),'FaceColor',kolory(:,j))
    plocik.FaceAlpha=0.9;
    hold on
    
    r0temp=cell2mat(r0(j));
    Rtemp(p)=Rr(j);
    G(r0temp==0)=NaN;
    AA(r0temp==0)=NaN; % rysowanie tylko tych przypadków, dla których istnieje orbita stacjonarna
    C=r0temp*100;
    Gz=G(I_Amin:I_Amax, I_gmin:I_gmax);
    Dz=(2*Con.nu./Gz).^(1/2)*skala;
    Cz=C(I_Amin:I_Amax, I_gmin:I_gmax);
    AAz=AA(I_Amin:I_Amax, I_gmin:I_gmax);
    
    %for s=1:3
    %contour(Gz,AAz,Cz,v1(s),'Color',kolory(:,j),'ShowText','off','Linewidth',3)%,'LineStyle',linest{s})
    [Ci,h]=contour(Dz,AAz,Cz,v1,'Color',kolory(:,j)/2,'ShowText','off','Linewidth',3,'LineStyle','--')
    clabel(Ci,h,'FontSize',13,'FontName','Helvetica','Margin',1,'LabelSpacing',280)
    hold on
    %end
    
end

%% Add simulation points plot
x=[1,1,0.2];
y=[0.0003,0.0014,0.003];
scatter(x,y,70,'filled','MarkerFaceColor','black')
tekst={'1)','2)','3)'};
for s=1:3
text(0.97*x(s),1.25*y(s),tekst(s),'FontSize',17,'FontWeight','bold','FontName','Helvetica')
hold on
end
%% Add line of A_cr
plot(fliplr(delta*skala),zeros(numel(delta),1)+Con.Acr,'k','LineWidth',2,'LineStyle','-.','Color',[0.2 0.2 0.2])
hold on

legend({'R=23um','','R=13um','','R=3um','','simulations','A_{cr}'},'FontSize',13,'FontName','Helvetica')
hold on

%% Set figure properties
xlim([delta_min*skala delta(1)*skala])
ylim([Amin 0.03])
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'fontsize',20)
xlabel('$\delta [cm]$','interpreter','latex','FontSize',30)
ylabel('$A$','interpreter','latex','FontSize',30)
set(gcf,'paperunits','centimeters')
set(gcf,'papersize',[width,height])
hold off
