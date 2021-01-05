clear
close all
clc
delete(gcp('nocreate'))
tic

%% Data
DIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Cloud_voids_theoretically/Stare_przydatne/';
fDIR=[DIR 'Functions/'];

% load spec_dir or spec_dir_even file
sciezki='spec_dir';
load([DIR,'Data/','Simulation_src_folders_ordered_20210104.mat'],sciezki)

if strcmp(sciezki,'spec_dir_even')
spec_dir=spec_dir_even;
end
% kolejnosc rysowania wg rozmiaru pliku trajektorii -numerki odnosza sie do
% nr-ow nadanych w tabelce excela
kolejka=[28];%,25,28];



%% Plot parameters
poolnr=12;
lcz_usr=50;
n_bin=75;
n_cz=5; % number of size parts for the plot
kat=40; % arbitrary Mie scattering angle
Rx_conf_value=0.95;
bright_Rx=1;
threshold=100/256;
fsize=13; % font size in the plot
width=16;
height=20;
skala=100; % to have a plot in [cm]
kolory=[0.4*linspace(2,5,n_cz)/n_cz; 0.8*linspace(2,5,n_cz)/n_cz; 1*linspace(2,5,n_cz)/n_cz]';
ymax=130*10^6; %50,130
%% Declare constants
Const = Constants;
%% Load data and functions
addpath(fDIR)
for jj=1:numel(kolejka)
    
    spec=spec_dir(kolejka(jj));
    loadDIR=char(strcat(DIR,spec_dir(kolejka(jj))));
    nazwa=char(strcat(DIR,spec_dir(kolejka(jj)),'histdane_lcz_usr_', num2str(lcz_usr), '_n_bin_', num2str(n_bin), '_n_cz_', num2str(n_cz), '.mat'));
    
    %% Calculation
    if isfile(nazwa)==0
        load([loadDIR,'Trajectories.mat'])
        load([loadDIR,'Parameters.mat'],'parameters')
        
        
        %% Times for averaging
        tplot_ind=numel(tglobal); % time or times at which the calculation is being made
        tglob_ind=1:numel(tglobal);
        for j=1:numel(tplot_ind)
            [~,I]=sort(abs(tglob_ind-tplot_ind(j)));
            czasy_ind=tglob_ind(I(1:lcz_usr));
            wynik{j}.czasy_ind=sort(czasy_ind);
        end
        
        % For some time R was not saved for polidispersed case
        if isfield(part(1).new_par,'R')==0
            for k=1:part(1).par.l_krop
                tau_p(k)=part(k).new_par.tau_p;
            end
            Rall=sqrt(9*Const.nu*Const.ro_a*tau_p/(2*Const.ro_p));
        else
            for k=1:part(1).par.l_krop
                Rall(k)=part(k).new_par.R;
            end
        end
        
        [Rx]=cut_dist(part(1).par.ProbDist,Rx_conf_value);
        
        [edges]=divide_dist(part(1).par.ProbDist,n_cz);
        
        [Iscat,Ix]=MieScatScaling(Rall*10^6,Rx*10^6,kat);
        col=Iscat*bright_Rx/Ix;
        col(col>1)=1; % opposite to what is drawn, col=1-col;
        
        for j=1:numel(tplot_ind)
            czasiki=wynik{j}.czasy_ind;
            polr=[];
            rozmR=[];
            bright=[];
            for l=1:numel(czasiki)
                t_ind=czasiki(l);
                droplets=drop_in_time{t_ind};
                r=zeros(1,numel(droplets));
                R=zeros(1,numel(droplets));
                Col=zeros(1,numel(droplets));
                for k=1:numel(droplets)
                    q=droplets(k);
                    a=find(abs(part(q).traj.t-tglobal(t_ind))<0.5*delta_t);
                    r(k)=sqrt(part(q).traj.X(a)^2+part(q).traj.Y(a)^2);
                    R(k)=Rall(q);
                    Col(k)=col(q);
                end
                polr=[polr,r];
                rozmR=[rozmR,R];
                bright=[bright,Col];
            end
            wynik{j}.polr=polr;
            wynik{j}.rozmR=rozmR;
            wynik{j}.bright=bright;
            
            [divis,binedges]=bin_division(n_bin,0,part(1).par.D,polr);
            bincent=binedges+0.5*(binedges(2)-binedges(1));
            bincent(end)=[];
            wynik{j}.radii_binedges=edges;
            wynik{j}.pos_bincenters=bincent;
            meanR=zeros(n_bin,1);
            meanRMie=zeros(n_bin,length(threshold));
            for m=1:n_bin
                meanR(m)=mean(rozmR(divis{m}));
                meanR(isnan(meanR)==1)=0;
                for p=1:numel(threshold)
                    meanRMie(m,p)=mean(nonzeros(rozmR(divis{m}).*(bright(divis{m})>threshold(p))));
                end
                meanRMie(isnan(meanRMie)==1)=0;
                wynik{j}.meanR=meanR;
                wynik{j}.meanRMie=meanRMie;
            end
            
            obj=2*pi*(binedges(2)-binedges(1))*(bincent)*2*part(1).par.Z;
            
            if sum(obj)-2*part(1).par.Z*pi*(part(1).par.D^2)>10^(-7)
                error('Cos nie tak z wyliczeniem objetosci')
            end
            
            for m=1:n_cz
                h=histogram(nonzeros((rozmR<=edges(m+1)).*polr),binedges);
                concentrations{m}=h.Values./(lcz_usr*obj);
                max_conc(m)=max(concentrations{m});
            end
            close all
            max_con(j)=max(max_conc);
            wynik{j}.n=concentrations;
            
        end
        ylim=max(max_con);
        
        %save(nazwa,'wynik')
        save(nazwa,'wynik','tplot_ind','n_cz','lcz_usr','ylim','-v7.3')
        clear part
    else
        load(nazwa)
        load([loadDIR,'Trajectories.mat'],'tglobal','drop_in_time')
        load([loadDIR,'Parameters.mat'],'parameters')
        tplot_ind=numel(tglobal);
    end
    
%% Plot figure
figure('Position', [600, 500, 1000, 8000])
set(gcf,'paperunits','centimeters')
set(gcf,'papersize',[width,height])


for j=1:numel(tplot_ind)
    for k=fliplr(1:n_cz)
        yyaxis left
        b=bar(wynik{j}.pos_bincenters*skala,wynik{j}.n{k}/10^6,'hist');
        b.FaceColor=kolory(k,:);
        hold on
    end
    set(gca,'xlim',[0 parameters.D*skala],...
        'ylim',[0 ymax/10^6],...
        'FontSize',fsize)
    xlabel('r [cm]','FontSize',fsize)
    xticks([0:1:parameters.D*skala])
    ylabel('n [1/cm^3]','FontSize',fsize)
    yticks(0:5:ymax/10^6)
    hold on
    
    yyaxis right
    plot(wynik{j}.pos_bincenters*skala, wynik{j}.meanR*10^(6),'Linewidth',3,'Color',[1 0.4 0.4])
    set(gca,'ylim',[0 15])
    yticks(0:2:15)
    hold on
    
    yyaxis right
    plot(wynik{j}.pos_bincenters*skala, wynik{j}.meanRMie(:,1)*10^(6),'Linewidth',3,'Color',[1 0.4 0.4])
    %set(gca,'ylim',[0 ylim])
    ylabel('<R> [\mu m], <S> [\mu m]','FontSize',fsize)
    hold off
    
    % plot(1:n_bin,meanR)
    % hold on
    % for p=1:numel(threshold)
    % plot(1:n_bin,meanRMie(:,p))
    % hold on
    % end
    % legend('Standard','1','2','3','4')
    
end
saveas(gcf,[loadDIR,'n_r_Rvis_mean_plot'],'fig')
saveas(gcf,[loadDIR,'n_r_Rvis_mean_plot'],'png')
fprintf(['Calculation for ' num2str(jj) ' finished.'])
end



function [division, binedges]=bin_division(binnr,edgeleft,edgeright,values)
%% Devides a set of values into binnr number of bins in the range [edgeleft,edgeright]
binedges=linspace(edgeleft,edgeright,binnr+1);
for m=1:binnr
    division{m}=find((values>binedges(m)).*(values<=binedges(m+1)));
end
end

function [edges]=divide_dist(PD,n)

D = cumsum(PD(:,2));
edge_val=(1:(n-1))/n;
edges=zeros(n,1);
for j=1:(n-1)
    [~,I]=min(abs(D-edge_val(j)));
    edges(j+1)=PD(I,1);
end
edges=[edges;1];

end

function [boarderline]=cut_dist(PD,conf_value)
D = cumsum(PD(:,2));
[~,I]=min(abs(D-conf_value));
boarderline=PD(I,1);
end