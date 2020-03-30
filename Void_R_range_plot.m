clear
close all
clc
delete(gcp('nocreate'))

DIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Cloud_voids/';

%% Load constants and functions
addpath(DIR)
Const=Constants;
fDIR=[DIR 'Functions/'];
addpath(fDIR)

%% Data
delta=[0.1:0.005:1.005]*10^(-2);
A=0.00001:2*10^(-5):0.03001;

%% Calulations

% R1

%R2 - A<A_cr

%R2 - A>A_cr

% Plot