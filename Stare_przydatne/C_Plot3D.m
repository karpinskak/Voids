%% Plot3D
% This script turnes the trajectory data around (if needed)
% and plots it in 3D with vectors showing particle velocity
% together with gravity vector.

clear
close all
clc
delete(gcp('nocreate'))
tic

%% Data
DIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Cloud_voids/';
spec_dir={
    'Results/A00014/Delta001/R12_9/Th05/k5_l2_n10_tk8_even_poli_exper/',...
    };

%% Declare constants
Const = Constants;

fDIR=[DIR 'Functions/'];
addpath(fDIR)

%% Parameters
% General
poolnr=1;
zapis=1; % 0 - do not save 1- save
% Droplet plot
co_ktora=10;
co_ile_s=0.005; % [s]
Mie=2; % particle size scaling 0 - no scaling, 1- linear, 2 - Mie,
ogon=1; % tail drawing: 0- no, 1- yes
czas=1; % 1- plot animation time, 0- do not plot animation time
if Mie==2
    kat=40;
    Rx_conf_value=0.95;
end

Rs=12.9*10^(-6);
size_Rs=7;
bright_Rx=1;

% Figure
tail_ratio=0.055; % czesc sekundy w czasie ktorej ruch rysujemy w postaci ogonka
campos_vec=[-0.64,-0.38, 0.11];
fsize=13; % font size in the plot
width=16;
height=20;
skala=100; % to have a plot in [cm]
dimension=3;


%% Plot
for jj=1:numel(spec_dir)
    %% Load data
    spec=spec_dir(jj);
    loadDIR=char(strcat(DIR, spec));
    load([loadDIR,'Trajectories.mat'],'tglobal')
    nazwa=[loadDIR,'Trajectories_rot.mat'];
    load([loadDIR,'Parameters.mat'],'dispersity','par_set','T')
    
    start=11551;
    stop=numel(tglobal);
    %start=1;
    %stop=6081;
    co_ile=round(co_ile_s*numel(tglobal)/T);
    
    %% Rotate
    if isfile(nazwa)==1
        load(nazwa)
        load([loadDIR,'Trajectories.mat'],'drop_in_time','delta_t')
    else
        load([loadDIR,'Trajectories.mat'])
        load([loadDIR,'Parameters.mat'],'T')
        % Parametry obrotu
        kat_obr=part(1).par.teta*360/(2*pi);
        kier_obr=[1,0,0];
        punkt_obr=[0,0,0];
        if poolnr>1
            parpool('local',poolnr)
        end
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
            colors=zeros(l_krop,1)+0.5;
            sizes=repmat(zeros(1,l_krop)+size,3);
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
    figure1 = figure('Color',[1 1 1]);
    Plot_features

    czasy_obr=start:co_ile:stop;
    
    if poolnr>1
        parpool('local',poolnr)
    end
    
    if czasy_obr(end)~=numel(tglobal)
        czasy_obr=[czasy_obr,numel(tglobal)];
    end
    
    for m=1:numel(czasy_obr)
        
        p=czasy_obr(m);
        
        droplets=cell2mat(drop_in_time(p));
        NumDrop=length(1:co_ktora:numel(droplets));
        kropelki=1:co_ktora:numel(droplets);
         Xrys=zeros(1,NumDrop);
         Yrys=zeros(1,NumDrop);
         Zrys=zeros(1,NumDrop);
         size_rys=zeros(1,NumDrop);
         colors_rys=zeros(NumDrop,3);
         Xtail=cell(1,NumDrop);
         Ytail=cell(1,NumDrop);
         Ztail=cell(1,NumDrop);
        for l=1:NumDrop
            ktora=kropelki(l);
            q=droplets(ktora);
            a=find(abs(traj(q).t-tglobal(p))<0.5*delta_t);
            Xrys(l)=skala*traj(q).X(a);
            Yrys(l)=skala*traj(q).Y(a);
            Zrys(l)=skala*traj(q).Z(a);
            size_rys(l)=sizes(q)^2;
            colors_rys(l,:)=colors(q,:);
            if a-kon(p)>=0
                Xtail{l}=skala*traj(q).X((a-kon(p)+1):a);
                Ytail{l}=skala*traj(q).Y((a-kon(p)+1):a);
                Ztail{l}=skala*traj(q).Z((a-kon(p)+1):a);
            else
                Xtail{l}=NaN;
                Ytail{l}=NaN;
                Ztail{l}=NaN;
            end
        end
        
        
        scatter3(Xrys', Yrys', Zrys', size_rys',colors_rys,'filled')
        hold on
        
        for ll=1:NumDrop
            plot3(Xtail{ll},Ytail{ll},Ztail{ll},'Color',[242/255,54/255,54/255])
            hold on
        end
        
       % clear size_rys colors_rys Xrys Yrys Zrys Xtail Ytail Ztail
        if czas==1
            czasik=0.1*floor(p*delta_t*10);
            str=['t=',num2str(czasik),'s'];
            text(0.9*xmax,0.9*ymin,0.9*zmax,str,'HorizontalAlignment','right','FontSize',fsize+3);
        end
        arrow3D([0,0,0], [0 0 -3.5],[ 0 0 1], 0.7)
        campos(campos_vec)
        Plot_features
        %hold on
        %     punkty_x=skala*punkty(:,1)*part(1).par.delta;
        %     punkty_y=skala*punkty(:,2)*part(1).par.delta;
        %     scatter(punkty_x,punkty_y,10);
        %     hold on
        
        %hold on
        [az,el] = view;
        view(az,el*1.0001)
        hold off
        if zapis==1
            nazwa=char(strcat(loadDIR, 'Animation_3D/Fig',num2str(p),'_3D'));
            saveas(gcf,nazwa,'fig')
            saveas(gcf,nazwa,'png')
        end
    end
end
delete(gcp('nocreate'))
function [boarderline]=cut_dist(PD,conf_value)
D = cumsum(PD(:,2));
[~,I]=min(abs(D-conf_value));
boarderline=PD(I,1);
end