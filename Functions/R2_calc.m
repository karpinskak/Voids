function [R_2] = R2_calc(Const,teta,A,delta,R_1,AA,Delta)
%R2_CALC Summary of this function goes here
%   Detailed explanation goes here
syms r
    assume(r,{'positive','real'})
    Sv_max=zeros(numel(A),1);
    for k=1:numel(A)
        A_temp=A(k);
        f_A=@(r) eq_curve(r,A_temp);
        if A_temp<Const.Acr
            df=diff(f_A,r);
            crit=double(vpasolve(df==0,r,[Const.rs,Const.ri]));
            max_id=find((crit<Const.ri).*(crit>Const.rs));
            if isempty(max_id)==0
                Svm=eq_curve(crit(max_id),A_temp);
            else
                error(['Rmax poza zakresem dla ',num2str(k)])
            end
        else
            Svm=eq_curve(Const.rs,A_temp);
        end
        Sv_max(k)=Svm;
    end
    SV_max=repmat(Sv_max,1,numel(delta));
    
    for l=1:numel(teta)
        tau_p=Const.nu*(Const.g*AA.*Delta*sin(teta(l))).^(-1).*SV_max;
        R2=10^6*sqrt(9*Const.ro_a*Const.nu*tau_p/(2*Const.ro_p));
        R1=R_1{l};
        R2(isnan(R1))=NaN;
        R_2{l}=R2;
    end
end

function krzywa=eq_curve(r,A)

chi=(1-exp(-r.^2/2))./(2*pi*A.*r.^2);
krzywa=A.*r.*sqrt(1+chi.^2);
end
