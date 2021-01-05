clear all
delta_tet05=[
0.01
0.005
0.008
0.004
0.008
0.006
0.002
0.002
];
A_tet05=[
0.0003
0.0003
0.0005
0.001
0.0018
0.002
0.003
0.006

];
delta_tet075=[
0.01
0.005
0.006
0.012
0.01
0.008
0.006
0.007

];
A_tet075=[
0.0003
0.0003
0.0008
0.0012
0.0014
0.0018
0.002
0.002

];
clf
figure(1)
scatter(delta_tet05*100,A_tet05)
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlim([0.1,1.2])
ylim([10^(-4),5*10^(-2)])
figure(2)
scatter(delta_tet075*100,A_tet075)
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlim([0.1,1.2])
ylim([10^(-4),5*10^(-2)])