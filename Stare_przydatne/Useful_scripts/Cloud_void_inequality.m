clear
clc
close all

%  syms M C D
%  assume(M > 0)
%  assume(C > 0)
% assume(D > 0)


%  C=[0:(10^(-6)):(2*10^(-5))].^2;
%   %C=10^(-5);
%   M=0:(10^(-6)):(4*10^(-4));
%   M=sort([-M,M]);
%   for j=1:numel(C)
%       
%   funk=M.^3-16*pi^2*C(j)*M-C(j);
%   plot(M, funk)
  %loglog(M,funk)
%   hold on
%   end

% %Amax=solve(M.^3-D*C*M-C,M,'MaxDegree',3)
%  %D=(4*pi)^2;
%  p=(C/2 + (- (D^3*C.^3)/27 + C.^2/4).^(1/2)).^(1/3) + (C*D)./(3*(C/2 + (C.^2/4 - (C.^3*D^3)/27).^(1/2)).^(1/3));
% % 
%  %plot(C,p)
% 
%  taylor(p,C) %'Order',5)%,'OrderMode', 'relative')


%%
% B=0:(10^(-7)):10^(-5);
% D=4*pi;
% A =sqrt((9*B.^2 + sqrt(3)*sqrt(27*B.^4 - 4*B.^6*D^6)).^(1/3)/(2^(1/3)*3^(2/3)) + ((2/3)^(1/3)*B.^2*D^2)./(9*B.^2 + sqrt(3)*sqrt(27*B.^4 - 4*B.^6*D^6)).^(1/3));
% Aap=2^(-1/6)*B.^(1/3);
%  plot(B,A,'r',B,Aap,'g')

% %%
%   syms B D
%   %D=4*pi;
%  A =sqrt((9*B^2 + sqrt(3)*sqrt(27*B^4 - 4*B^6*D^6))^(1/3)/(2^(1/3)*3^(2/3)) + ((2/3)^(1/3)*B^2*D^2)/(9*B^2 + sqrt(3)*sqrt(27*B^4 - 4*B^6*D^6))^(1/3));
%  xxx=simplify(series(A,B,0,'Order',1));

%%
r=0.001:0.1:1.5;
AA=[0.0001,0.003,0.006,0.015,0.02176,0.03];
for j=1:numel(AA)
    A=AA(j);
krzywa_org=r*A.*(1+((1-exp(-r.^2/2))./(2*pi*A*r.^2)).^2).^(1/2);
krzywa_aprox=r*sqrt(1+(4*pi*A)^2)/(4*pi);%-r.^3*(16*pi*sqrt(1+(4*pi*A)^2))^(-1);
%krzywa_aprox=r/(4*pi);
plot(r,krzywa_org,'r')
hold on
plot(r,krzywa_aprox,'g')
hold off
pause(1)
end


% syms r A
%  taylor(r*A*(1+((1-exp(-r^2/2))/(2*pi*A*r^2))^2)^(1/2),r,0,'Order',4)