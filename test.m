
Const=Constants;
teta_temp=0.5*pi/2;
delta_temp=0.005;
syms tau_p 
assume(tau_p,{'real','positive'})
A_temp=0.015;
P=@(tau_p) 1+Const.nu^(-1)*Const.g^2*tau_p.^3.*(sin(teta_temp)).^2-Const.nu*tau_p.*delta_temp.^(-2);
Amax=@(tau_p) (0.5*(4*pi)^(-2)*P(tau_p).*((1+4*Const.nu*tau_p./(P(tau_p).^2.*delta_temp.^2)).^(1/2)-1)).^(1/2);
tau_p_sol=double(vpasolve(Amax(tau_p)-A_temp==0,tau_p));
R_sol=10^6*sqrt(9*Const.ro_a*Const.nu*tau_p_sol/(2*Const.ro_p));
abs(Amax(tau_p_sol)-A_temp)<10^(-7)
