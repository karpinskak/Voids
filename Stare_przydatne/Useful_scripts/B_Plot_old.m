clear
close all
clc
delete(gcp('nocreate'))
tic

%% Data
DIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Cloud_voids/';

spec_dir={...
    'Results/A0002/Delta0007/R12_9/Th05/k7_14_l2_n10_tk2_poli_exper/',...
    };
%025 dla Th oznacza 0.25*pi/2

fDIR=[DIR 'Functions/'];
%% Declare constants
Const = Constants;
%% Load data and functions
addpath(fDIR)

%% Animation parameters
dimension=2; % 2 or 3
anim=1; % 0 - for whole trajectories, 1 - on
zapis=0; % 0 - do not save 1- save
co_ktora=500;
punkty_rys=1; % plot equillibrium points 1- yes, 0 - no
if anim==1
    co_ile=1000;
    start_d='end'; % somehow general times
    stop_d='end'; %
    Mie=2; % particle scaling 0 - no scaling, 1- brightness linear, 2 - Mie,
    czas=1; % 0- plot animation time, 1- do not plot animation time
    sl_width=4/100; % [m] plot particles from the central slice of the domain [z_middle-sl_width,z_middle+sl_width]
    z_middle=1/100;
    if Mie==2
        kat=40;
        Rx_conf_value=0.95;
    end
end

%% Plot parameters
Rs=12.9*10^(-6);
size_Rs=40;
bright_Rx=1;
fsize=13; % font size in the plot
width=16;
height=20;
skala=100; % to have a plot in [cm]

%% Plot
%for jj=1:numel(spec_dir)
    jj=1;
    % Load data
    spec=spec_dir(jj);
    loadDIR=char(strcat(DIR,spec_dir(jj)));
    load([loadDIR,'Trajectories.mat'])
    load([loadDIR,'Parameters.mat'],'dispersity','par_set')
    start=numel(tglobal);
    stop=numel(tglobal);
    
    %% Figure ranges
    xmin=-1.1*part(1).par.D*skala;
    xmax=1.1*part(1).par.D*skala;
    ymin=-1.1*part(1).par.D*skala;
    ymax=1.1*part(1).par.D*skala;
    zmin=-1.1*part(1).par.Z*skala;
    zmax=1.1*part(1).par.Z*skala;
    
    %% calculate size and color
    
    
    if anim==1
        if dispersity~=0
            [Rx]=cut_dist(part(1).par.ProbDist,Rx_conf_value);
        end
        
        switch Mie
            case 0
                size=part(1).par.R*size_Rs/Rs;
                colors=[zeros(part(1).par.l_krop,1)+0.5 zeros(part(1).par.l_krop,1)+0.5 zeros(part(1).par.l_krop,1)+0.5];
                sizes=zeros(1,part(1).par.l_krop)+size;
            case 1
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
                sizes=Rall*size_Rs/Rs;
                col=Rall*bright_Rx/Rx;
                col(col>=1)=1;
                col=1-col;
                colors=repmat(col,3);
            case 2
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
    end
    
    %% Calculate equillibrium points positions
    if dispersity==2||dispersity==3
        [bl]=cut_dist(part(1).par.ProbDist,0.95);
        promien=10^(-6):10^(-6):bl;
        punkty=[];
        stabilnosc=[];
        for k=1:numel(promien)
            [par_eq]=wylicz_param(Const,par_set,promien(k),part(1).par.delta,part(1).par.teta,part(1).par.A,1,1);
            [points, stability] = eq_points(part(1).par.A,par_eq.Sv,par_eq.St);
            punkty=[punkty;points];
            stabilnosc=[stabilnosc,stability];
        end
    end
    
%% Create figure
    figure1 = figure('Color',[1 1 1]);
    Plot_features
    
    if anim==1
        czasy_obr=start:co_ile:stop;
        for m=1:numel(czasy_obr)
            
            p=czasy_obr(m);
            
            droplets=cell2mat(drop_in_time(p));
            
            for l=1:co_ktora:numel(droplets)
                
                q=droplets(l);
                a=find(abs(part(q).traj.t-tglobal(p))<0.5*delta_t);
                if abs(part(q).traj.Z-z_middle)<=sl_width
                    if dimension==3
                        scatter3(skala*part(q).traj.X(a), skala*part(q).traj.Y(a), skala*part(q).traj.Z(a),sizes(q),'MarkerEdgeColor',colors(q,:),'MarkerFaceColor',colors(q,:))
                    elseif dimension==2
                        scatter(skala*part(q).traj.X(a), skala*part(q).traj.Y(a),sizes(q),'MarkerEdgeColor',colors(q,:),'MarkerFaceColor',colors(q,:));
                    end
                end
                hold on
                Plot_features
            end
            if czas==1
                czas=0.1*floor(p*delta_t*10);
                str=['t=',num2str(czas),'s'];
                if dimension==3
                    text(skala*part(1).par.D,skala*part(1).par.D,skala*0.09*part(1).par.Z,str,'HorizontalAlignment','right','FontSize',fsize+3);
                elseif dimension==2
                    text(skala*part(1).par.D,skala*part(1).par.D,str,'HorizontalAlignment','right','FontSize',fsize+3);
                end
            end
            hold on
            if punkty_rys==1
            punkty_x=skala*punkty(:,1)*part(1).par.delta;
            punkty_y=skala*punkty(:,2)*part(1).par.delta;
            if dimension==2
                scatter(punkty_x,punkty_y,10);
                hold on
            end
            end
            
            pause(1)
            hold off
            if zapis==1
                nazwa=char(strcat(loadDIR, num2str(p)));
                saveas(gcf,nazwa,'fig')
            end
        end
    elseif anim==0
        for j=1:co_ktora:l_krop
            xplot=part(j).traj.X;
            yplot=part(j).traj.Y;
            if isempty(xplot)==0 && isempty(yplot)==0
                kolor=[1; 0; 0]*mod(abs(part(1).par.D-part(j).traj.X(end))/(2*part(1).par.D),1)+...
                    [0 ;1 ;0]*mod(abs(-part(1).par.D-part(j).traj.X(end))/(2*part(1).par.D),1)+...
                    [0; 0; 1]*mod(abs(part(1).par.D-part(j).traj.Y(end))/(2*part(1).par.D),1)+...
                    [0; 0; 0]*mod(abs(-part(1).par.D-part(j).traj.Y(end))/(2*part(1).par.D),1);
                
                testowa=((xplot==0).*(yplot==0));
                xplot(testowa==1)=[];
                yplot(testowa==1)=[];
                xplot=[part(j).init.r0.*cos(part(j).init.fi0);xplot];
                yplot=[part(j).init.r0.*sin(part(j).init.fi0);yplot];
                plot(xplot,yplot,'Color',kolor)
            end
            hold on
        end
        punkty_x=punkty(:,1)*part(1).par.delta;
        punkty_y=punkty(:,2)*part(1).par.delta;
        scatter(punkty_x,punkty_y,10);
        if zapis==1
            nazwa=char(strcat(loadDIR, 'trajectories_all'));
            saveas(gcf,nazwa,'fig')
        end
        
        hold off
    end
    
%end

function [boarderline]=cut_dist(PD,conf_value)
D = cumsum(PD(:,2));
[~,I]=min(abs(D-conf_value));
boarderline=PD(I,1);
end
