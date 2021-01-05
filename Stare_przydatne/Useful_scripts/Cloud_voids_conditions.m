clear
clc

%% Declare constants
Const = Constants;

[delta,teta]=meshgrid(0.001:0.0001:0.02,pi/8:pi/16:pi/2);

B=(Const.rs*Const.nu^2/(Const.g*2^8*pi^3))*((sin(teta).*delta.^3).^(-1));
Amax=2^(-1/6)*B.^(1/3);
St_max_min=2^4*pi^2*Amax;

tau_p=delta.^2.*St_max_min.*Amax/(Const.nu);
R=10^6*(9*Const.nu*tau_p*Const.ro_a/(2*Const.ro_p)).^(1/2);

gamma=2*Const.nu*delta.^(-2);
figure(1)
subplot(2,1,1)
mesh(delta*100,teta*360/(2*pi),Amax)
xlabel('\delta [cm]')
ylabel('\theta [deg]')
zlabel('A_{max}')
subplot(2,1,2)
mesh(delta*100,teta*360/(2*pi),R)
xlabel('\delta [cm]')
ylabel('\theta [deg]')
zlabel('max(R_{min}')

figure(2)
loglog(gamma(1,:),Amax)
xlim([10^(-1) 10^(1)])
ylim([10^(-6) 10^(-2)])
xlabel('\gamma [1/s]')
ylabel('A_{max}')