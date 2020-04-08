function [R_m] =Rm_calc(Const,teta,A,delta)

for l=1:numel(teta)
    teta_temp=teta(l);
    Rm=zeros(numel(delta),numel(A));
    
    parfor j=1:numel(delta)
        delta_temp=delta(j);
        Rm_slice=zeros(1,numel(A));
        tau_p=sym('tau_p');
        P=@(tau_p) 1+Const.nu^(-1)*Const.g^2*tau_p.^3.*(sin(teta_temp)).^2-Const.nu*tau_p.*delta_temp.^(-2);
        Amax=@(tau_p) (0.5*(4*pi)^(-2)*P(tau_p).*((1+4*Const.nu*tau_p./(P(tau_p).^2.*delta_temp.^2)).^(1/2)-1)).^(1/2);
        for k=1:numel(A)
            A_temp=A(k);
            tau_p_sola1=double(vpasolve(Amax(tau_p)-A_temp==0,tau_p,5));
            tau_p_sola2=double(vpasolve(Amax(tau_p)-A_temp==0,tau_p,10));
            tau_p_sol=min(tau_p_sola1(tau_p_sola1>0),tau_p_sola2(tau_p_sola2>0));
           if isempty(tau_p_sol)==1 || tau_p_sol<2*Const.ro_p*(17*(10^(-6)).^2/(9*Const.nu*Const.ro_a))
               
               tau_p_sol=double(vpasolve(Amax(tau_p)-A_temp==0,tau_p,100));
           end
            R_sol=10^6*sqrt(9*Const.ro_a*Const.nu*tau_p_sol/(2*Const.ro_p));
            if isempty(tau_p_sol)==1 || imag(R_sol)~=0
                R_sol=NaN;
            else
            end
            Rm_slice(k)=real(R_sol);
        end
        Rm(j,:)=Rm_slice;
        disp([num2str(l),', ',num2str(j)])
    end
    R_m{l}=Rm';
    toc
end

end