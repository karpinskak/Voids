clear
close all
clc
delete(gcp('nocreate'))


%% Declare constants
Const = Constants;

%% Load functions
fDIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Cloud_voids_theoretically/Stare_przydatne/Functions/';
addpath(fDIR)

Sim_data_sets=table2struct(readtable('Dane_3.xls'));

for set_nr=1:numel(Sim_data_sets)
    tic
    [DIR,spec,poolnr,tstart,T,par_set,a,b,c,d,dispersity,Rmin,Rdev,PD,type,k,l,n,max_R,new_sim,loadPosDIR,loadDIR] = Parameters_sets(set_nr,Sim_data_sets);
    
    %% Create directory if not existing
    if isfolder(loadDIR)~=1
        mkdir(loadDIR)
    end
    
    %% Calculate all parameters needed
    [parameters]=wylicz_param(Const, par_set,a,b,c,d,k,l);
    clear a b c d k l
    part.par=parameters;
    switch dispersity
        case 2
            part.par.Rdev=Rdev;
            part.par.ProbDist=PD;
        case 1
            part.par.Rmax=Rmax;
        case 3
            part.par.ProbDist=PD;
    end
    part.par.Rmin=Rmin;
    
    save([loadDIR '/Parameters.mat'],'max_R','n','T','tstart','parameters','dispersity','par_set','type','new_sim')
    clear parameters
    fprintf('Parameters generated')
    
    %% Create initial conditions
    Generate_initial_conditions
    fprintf('Initial conditions generated')
    toc
    
    %% Calculate trajectories
    Trajectories
    fprintf(['Calculation for ' num2str(set_nr) ' finished.'])
    toc
    delete(gcp('nocreate'))
    clearvars -except Const Sim_data_sets
end

