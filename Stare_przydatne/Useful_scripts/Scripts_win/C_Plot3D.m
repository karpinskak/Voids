%% Plot3D
% This script turnes the trajectory data around (if needed) 
% and plots it in 3D with vectors showing particle velocity
% together with gravity vector.

clear
close all
clc
%delete(gcp('nocreate'))
tic

%% Data
DIR='C:\Users\ImaCo\Desktop\Kasia\Symulacje\';

spec={...
    'Results/A00014/Delta001/R12_9/Th05/k5_l2_n10_tk8_poli_exper/',...
    %'Results/A00003/Delta001/R12_9/Th075/k5_l2_n10_tk8_poli_exper/',...
    %'Results\A0003\Delta05\R12_9\Th075\k10_l2_n10_tk6_poli_exper\',...
    };

%% Declare constants
Const = Constants;
%% Load data
loadDIR=char(strcat(DIR, spec));
load([loadDIR,'Trajectories.mat'],'tglobal')

%% Parameters
% General
poolnr=3;
zapis=1; % 0 - do not save 1- save
% Droplet plot
co_ktora=10;
start=1;
stop=numel(tglobal);
 co_ile_s=0.01; % [s]
Mie=2; % particle size scaling 0 - no scaling, 1- linear, 2 - Mie,
ogon=1; % tail drawing: 0- no, 1- yes
czas=0; % 0- plot animation time, 1- do not plot animation time
if Mie==2
    kat=40;
    Rx_conf_value=0.95;
end

Rs=12.9*10^(-6);
size_Rs=40;
bright_Rx=1;

% Figure
tail_ratio=0.055; % czesc sekundy w czasie ktorej ruch rysujemy w postaci ogonka
campos_vec=[-0.64,-0.38, 0.11];
fsize=13; % font size in the plot
width=16;
height=20;
skala=100; % to have a plot in [cm]

%% Load data and functions
fDIR=[DIR 'Functions/'];
addpath(fDIR)
nazwa=[loadDIR,'Trajectories_rot.mat'];
load([loadDIR,'Parameters.mat'],'dispersity','par_set')
dimension=3;

%% Rotate
if  exist(nazwa, 'file') == 2
    load(nazwa)
    load([loadDIR,'Trajectories.mat'],'drop_in_time','delta_t')
else
    load([loadDIR,'Trajectories.mat'])
    load([loadDIR,'Parameters.mat'],'T')
    % Parametry obrotu
    kat_obr=part(1).par.teta*360/(2*pi);
    kier_obr=[1,0,0];
    punkt_obr=[0,0,0];
%     if poolnr>1
%         parpool('local',poolnr)
%     end
    parfor p=1:part(1).par.l_krop
        XYZold=[part(p).traj.X,part(p).traj.Y,part(p).traj.Z]';
        if isempty(part(p).traj.X)~=1
            [XYZnew, ~, ~] = AxelRot(XYZold, kat_obr, kier_obr, punkt_obr);
            traj(p).X=XYZnew(1,:);
            traj(p).Y=XYZnew(2,:);
            traj(p).Z=XYZnew(3,:);
        else
            traj(p).X=[];
            traj(p).Y=[];
            traj(p).Z=[];
        end
        traj(p).t=part(p).traj.t;
    end
%    delete(gcp('nocreate'))
    rot_parameters.angle=kat_obr;
    rot_parameters.direct=kier_obr;
    rot_parameters.point=punkt_obr;
    
    if dispersity~=0
        [Rx]=cut_dist(part(1).par.ProbDist,Rx_conf_value);
    end
    D=part(1).par.D;
    Z=part(1).par.Z;
    save(nazwa,'traj','rot_parameters','Rx','T','D','Z','-v7.3')
    switch Mie
        case 0
            R=part(1).par.R;
            l_krop=part(1).par.l_krop;
            save(nazwa,'R','l_krop','-append')
        otherwise
            if isfield(part(1).new_par,'R')==0
                for j=1:part(1).par.l_krop
                    tau_p(j)=part(j).new_par.tau_p;
                end
                Rall=sqrt(9*Const.nu*Const.ro_a*tau_p/(2*Const.ro_p));
            else
                for j=1:part(1).par.l_krop
                    Rall(j)=part(j).new_par.R;
                end
            end
            save(nazwa,'Rall','-append')
    end
end
clear part

%% calculate size and color

switch Mie
    case 0
        size=R*size_Rs/Rs;
        colors=[zeros(l_krop,1)+0.5 zeros(l_krop,1)+0.5 zeros(l_krop,1)+0.5];
        sizes=zeros(1,l_krop)+size;
    case 1
        sizes=Rall*size_Rs/Rs;
        col=Rall*bright_Rx/Rx;
        col(col>=1)=1;
        col=1-col;
        colors=repmat(col,3);
    case 2
        [Iscat,Ix]=MieScatScaling(Rall*10^6,Rx*10^6,kat);
        sizes=Rall*size_Rs/Rs;
        col=Iscat*bright_Rx/Ix;
        col(col>1)=1;
        col=1-col;
        if iscolumn(col)==0
            col=col';
        end
        colors=repmat(col,3);
end

%% Figure ranges
xmin=-1.7*D*skala;
xmax=1.7*D*skala;
ymin=-1.7*D*skala;
ymax=1.7*D*skala;
zmin=-1.1*Z*skala;
zmax=1.1*Z*skala;

%% Calculate tail lengths
%%% The tail is not drawn until the particles are not present long enough
dlug=floor(numel(tglobal)*tail_ratio/T); % length of the tail in the units of global time instants
kon=zeros(numel(tglobal),1);
kon((dlug+1):end,1)=dlug;

%% Create figure
set(gcf,'Renderer','OpenGL');
figure1 = figure('Color',[1 1 1]);
Plot_features

co_ile=round(co_ile_s*numel(tglobal)/T);
czasy_obr=start:co_ile:stop;
if czasy_obr(end)~=numel(tglobal)
    czasy_obr=[czasy_obr,numel(tglobal)];
end

for m=1:numel(czasy_obr)
    
    p=czasy_obr(m);
    
    droplets=cell2mat(drop_in_time(p));
    
    for l=1:co_ktora:numel(droplets)
        
        q=droplets(l);
        a=find(abs(traj(q).t-tglobal(p))<0.5*delta_t);
        if abs(skala*traj(q).Z(a))<=zmax && abs(skala*traj(q).X(a))<=xmax && abs(skala*traj(q).Y(a))<=ymax
            scatter3(skala*traj(q).X(a), skala*traj(q).Y(a), skala*traj(q).Z(a),sizes(q),'MarkerEdgeColor',colors(q,:),'MarkerFaceColor',colors(q,:))
            hold on
            if ogon==1
                if a-kon(p)>=0
                    Xtail=skala*traj(q).X((a-kon(p)+1):a);
                    Ytail=skala*traj(q).Y((a-kon(p)+1):a);
                    Ztail=skala*traj(q).Z((a-kon(p)+1):a);
                    plot3(Xtail,Ytail,Ztail,'Color',[242/255,54/255,54/255])
                end
            end
            hold on
        end
    end
    
    
    
    if czas==1
        czas=0.1*floor(p*delta_t*10);
        str=['t=',num2str(czas),'s'];
        text(skala*D,skala*D,skala*0.09*Z,str,'HorizontalAlignment','right','FontSize',fsize+3);
    end
    
    arrow3D([0,0,0], [0 0 -3.5],[ 0 0 1], 0.7)
    campos(campos_vec)
    Plot_features
    %hold on
    %     punkty_x=skala*punkty(:,1)*part(1).par.delta;
    %     punkty_y=skala*punkty(:,2)*part(1).par.delta;
    %     scatter(punkty_x,punkty_y,10);
    %     hold on
    %pause(1)
    hold off
    [az,el] = view;
    view(az,el*1.0001)
    if zapis==1
        nazwa=char(strcat(loadDIR, 'Animation_3D/Fig',num2str(p),'_3D'));
        saveas(gcf,nazwa,'fig')
        saveas(gcf,nazwa,'png')
    end
end
