
Const=Constants;
teta=pi/16;
delta=0.0011:0.00005:0.0015;
A=[0.0011,0.0123,0.0125,0.0127,0.0129];

[R_1] =R1_calc(Const,teta,A,delta);
 
function [R_1] =R1_calc(Const,teta,A,delta)

for l=1:numel(teta)
    teta_temp=teta(l);
    R1=zeros(numel(delta),numel(A));
    
    parfor j=1:numel(delta)
        delta_temp=delta(j);
        R1_slice=zeros(1,numel(A));
        tau_p=sym('tau_p');
        P=@(tau_p) 1+Const.nu^(-1)*Const.g^2*tau_p.^3.*(sin(teta_temp)).^2-Const.nu*tau_p.*delta_temp.^(-2);
        Amax=@(tau_p) (0.5*(4*pi)^(-2)*P(tau_p).*((1+4*Const.nu*tau_p./(P(tau_p).^2.*delta_temp.^2)).^(1/2)-1)).^(1/2);
        for k=1:numel(A)
            A_temp=A(k);
            tau_p_sol=double(vpasolve(Amax(tau_p)-A_temp==0,tau_p,[0,Inf]));
            R_sol=10^6*sqrt(9*Const.ro_a*Const.nu*tau_p_sol/(2*Const.ro_p));
            if isempty(tau_p_sol)==1 || imag(R_sol)~=0
                R_sol=NaN;
            else
            end
            R1_slice(k)=real(R_sol);
        end
        R1(j,:)=R1_slice;
        disp([num2str(l),', ',num2str(j)])
    end
    R_1{l}=R1';
end
end