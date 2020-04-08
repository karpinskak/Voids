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
nazwar1='R1_void_minimum';
nazwar2='R2_void_maximum';
nazwar3='Rm_void_middle';

%% Data
npool=13;
teta=[pi/16,pi/8,pi/4];
delta=[0.1:0.005:1.005]*10^(-2);
A=0.0001:2*10^(-4):0.03001;
%delta=[0.1:0.005:0.2]*10^(-2);
%A=0.0001;%:(2*10^(-4)):0.002;
[Delta,AA]=meshgrid(delta,A);


% Declarations
R=zeros(numel(delta),numel(A));
R_1=cell(numel(teta),1);
R_2=cell(numel(teta),1);

% Calculations
%% R1
if isfile([DIR nazwar1, '.mat'])
    load([DIR nazwar1,'.mat'])
else
    parpool('local',npool)
    tic
    [R_1] =R1_calc(Const,teta,A,delta);
    save([DIR nazwar1,'.mat'],'teta', 'delta', 'A', 'R_1', 'DIR')
    toc
end

%% Rm
if isfile([DIR nazwar3, '.mat'])
    load([DIR nazwar3,'.mat'])
else
    disp('Calculating Rm now.')5c
    parpool('local',npool)
    tic
    [R_m] =Rm_calc(Const,teta,A,delta);
    save([DIR nazwar3,'.mat'],'teta', 'delta', 'A', 'R_m', 'DIR')
    toc
end



%% R2 - A<A_cr
% finding Svmax first for given A in the range Svs Svi, then calculating R out of it.
%R2 - A>A_cr
% equillibrium curve value for rs and given A gives Sv, then calc R

if isfile([DIR nazwar2 '.mat'])
    load([DIR nazwar2,'.mat'])
else
    [R_2] = R2_calc(Const,teta,A,delta,R_1,AA,Delta);
end
delete(gcp('nocreate'))
save([nazwar2,'.mat'],'teta', 'delta', 'A', 'R_2', 'DIR')

%% Plot

% Full plot for thesis
Void_R_range_full_plot


% % test plot
% tau_p=10^(-7):10^(-3):5;
% teta_temp=teta;
% for j=1:numel(A)
%     A_temp=A(j);
%      for k=1:numel(delta)
%          delta_temp=delta(k);
%         P= 1+Const.nu^(-1)*Const.g^2*tau_p.^3.*(sin(teta_temp)).^2-Const.nu*tau_p.*delta_temp.^(-2);
%         Amax_minA=(0.5*(4*pi)^(-2)*P.*((1+4*Const.nu*tau_p./(P.^2.*delta_temp.^2)).^(1/2)-1)).^(1/2)-A_temp;
%     plot(tau_p,Amax_minA)
%     hold on
%      end
% end