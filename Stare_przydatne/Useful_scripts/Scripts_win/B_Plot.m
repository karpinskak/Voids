clear
close all
clc
%delete(gcp('nocreate'))
tic

%% Data
DIR='C:\Users\ImaCo\Desktop\Kasia\Symulacje\';

spec_dir={...
    %'Results/A00014/Delta001/R12_9/Th05/k5_l2_n10_tk8_poli_exper/',...
    %'Results/A00003/Delta001/R12_9/Th075/k5_l2_n10_tk8_poli_exper/',...
    'Results/A0003/Delta05/R12_9/Th075/k10_l2_n10_tk6_poli_exper/',...
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
zapis=1; % 0 - do not save 1- save
co_ktora=1;
punkty_rys=0; % plot equillibrium points 1- yes, 0 - no
if anim==1
    co_ile_s=0.01; % [s]
    start_d=1; % somehow general times
    stop_d='end'; %
    Mie=2; % particle scaling 0 - no scaling, 1- brightness linear, 2 - Mie,
    czas=1; % 0- plot animation time, 1- do not plot animation time
    sl_width=4/100; % [m] plot particles from the central slice of the domain [z_middle-sl_width,z_middle+sl_width]
    z_middle=0/100;
    if Mie==2
        kat=40;
        Rx_conf_value=0.95;
    end
end

%% Plot parameters
Rs=12.9*10^(-6);
size_Rs=0.07; % [cm]
bright_Rx=1;
fsize=13; % font size in the plot
width=16;
height=20;
skala=100; % to have a plot in [cm]

%% Plot
for jj=1:numel(spec_dir)

    % Load data
    spec=spec_dir(jj);
    loadDIR=char(strcat(DIR,spec_dir(jj)));
    load([loadDIR,'Trajectories.mat'])
    load([loadDIR,'Parameters.mat'],'dispersity','par_set','T')
    start=1;
    stop=numel(tglobal);
    co_ile=round(co_ile_s*numel(tglobal)/T);
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
                %colors=repmat(col,3);
                colors=col;
        end
    end
    
    %% Calculate equillibrium points positions
    if punkty_rys==1
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
    end
%% Create figure
    figure1 = figure('Color',[1 1 1]);
    Plot_features
    
    if anim==1
        czasy_obr=start:co_ile:stop;
        if czasy_obr(end)~=numel(tglobal)
            czasy_obr=[czasy_obr,numel(tglobal)];
        end
            
        for m=1:numel(czasy_obr)
            
            p=czasy_obr(m);
            
            droplets=cell2mat(drop_in_time(p));
      
            licz=0;
            
            for l=1:co_ktora:numel(droplets)
                
                q=droplets(l);
                a=find(abs(part(q).traj.t-tglobal(p))<0.5*delta_t);
               
                
                if abs(part(q).traj.Z-z_middle)<=sl_width
                    licz=licz+1;
                    Xrys(licz)=skala*part(q).traj.X(a);
                    Yrys(licz)=skala*part(q).traj.Y(a);
                    size_rys(licz)=sizes(q);
                    colors_rys(licz)=colors(q);
                    if dimension==3
                        Zrys(licz)=skala*part(q).traj.Z(a);
                    end
                end
            end
            NumDrop=numel(size_rys);
            pkt = linspace(0,2*pi,32)';
            pkt_okregu = [cos(pkt) sin(pkt)];
            colormap('gray')
            patch(bsxfun(@plus,bsxfun(@times,size_rys,repmat(pkt_okregu(:,1),[1 NumDrop])),Xrys),...
                bsxfun(@plus,bsxfun(@times,size_rys,repmat(pkt_okregu(:,2),[1 NumDrop])),Yrys),...
                repmat(colors_rys,32,1),...
                'EdgeColor','none');
                
                hold on
                Plot_features
                if punkty_rys==1
                    punkty_x=skala*punkty(:,1)*part(1).par.delta;
                    punkty_y=skala*punkty(:,2)*part(1).par.delta;
                    if dimension==2
                        scatter(punkty_x,punkty_y,10);
                        hold on
                    end
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
               
                
                hold off
                if zapis==1
                    nazwa=char(strcat(loadDIR, 'Fig', num2str(p), '_2D'));
                    saveas(gcf,nazwa,'fig')
                    saveas(gcf,nazwa,'png')
                end
                clf
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
            saveas(gcf,nazwa,'png')
        end
        
        hold off
        clf
    end
    
end