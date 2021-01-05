function [R_s] = Rs_calc(Const,teta,A,delta,AA,Delta)
%R2_CALC Summary of this function goes here
%   Detailed explanation goes here

    Sv_max=zeros(numel(A),1);
    for k=1:numel(A)
        A_temp=A(k);
        Svm=eq_curve(Const.rs,A_temp);
        Sv_max(k)=Svm;
    end
    SV_max=repmat(Sv_max,1,numel(delta));
    
    for l=1:numel(teta)
        tau_p=Const.nu*(Const.g*AA.*Delta*sin(teta(l))).^(-1).*SV_max;
        Rs=10^6*sqrt(9*Const.ro_a*Const.nu*tau_p/(2*Const.ro_p));
        %R1=R_1{l};
        %R2(isnan(R1))=NaN;
        R_s{l}=Rs;
    end
end

function krzywa=eq_curve(r,A)

chi=(1-exp(-r.^2/2))./(2*pi*A.*r.^2);
krzywa=A.*r.*sqrt(1+chi.^2);
end
