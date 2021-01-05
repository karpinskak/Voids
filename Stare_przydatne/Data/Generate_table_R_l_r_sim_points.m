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

nazwar1='R1_void_sim_points';
nazwar2='R2_void_sim_points';
nazwar3='Rs_void_sim_points';

if isfile([DIR nazwar1, '.mat']) && isfile([DIR nazwar2, '.mat']) && isfile([DIR nazwar3, '.mat'])
    load([DIR nazwar1,'.mat'])
    load([DIR nazwar2,'.mat'])
    load([DIR nazwar3,'.mat'])
else
    error('Lacking data file.')
end
load([DIR,'Stare_przydatne/Data/Simulation_vortex_parameters.mat'])


%% Calculate R_<, R_r, deltaR

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

%% Find simulation points in grid indices
for jj=1:numel(vortex_param.A)
j=find(vortex_param.A(jj)==A);
k=find(vortex_param.delta(jj)==delta);
l=find(vortex_param.teta(jj)==teta);
R_l_r_sim_points(jj,1)=R_l{l}(j,k);
R_l_r_sim_points(jj,2)=R_r{l}(j,k);

end

save([DIR,'Stare_przydatne/Data/Simulation_vortex_parameters.mat'],'vortex_param','R_l_r_sim_points','nazwar1','nazwar2','nazwar3')
