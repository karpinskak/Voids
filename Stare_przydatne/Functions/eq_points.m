function [points, stability] = eq_points(A,Sv,St)
%EQ_POINTS Calculates all eq. points of Burgers vortex and its stability
%   'points' is nx2 matrix where n is the number of eq. points and 2 gives X
%   and Y positions.
%   'stability' can have values 1,2 which are stable and unstable
% 20210106: 1) Part for A>Acr was wrong - "points" had just
% radial component calculated, stability was assumed =1. I  just finished
% the if clause earlier and applied the calculations of x,y point and it
% stability as for the other case.

syms r
assume(r>=0)

if A>=0.02176
    point_r=vpasolve(r*sqrt(1+((1-exp(-r^2/2))/(2*pi*A*r^2))^2)-Sv/A==0,r);

else
    rmax=200;
    r_bezw=[0:0.001:6,6.01:0.01:rmax];
    [~,kolejnosc]=sort(abs(krzywa(r_bezw,A,Sv)));
    r_pos=r_bezw(kolejnosc);
    r0=[r_pos(1:10),rmax,10*rmax];
    for n=1:numel(r0)
        x0=r0(n);
        r_equ=double(vpasolve(r.*sqrt(1+((1-exp(-r.^2/2))./(2*pi*A*r.^2)).^2)-Sv/A==0,r,x0));
        if numel(r_equ)==1
            r_eq0(n)=r_equ;
        elseif numel(r_equ)==0
            r_eq0(n)=NaN;
        else
            error('za duzo rozwiazan')
        end
    end
    r_eq0(isnan(r_eq0))=[];
    point_r=unique(r_eq0);
    if numel(point_r)>3
        error(['Sth wrong with nr of eqillibrium points:' num2str(numel(points))])
    end
end

    points=zeros(numel(point_r),2);
    chi=(1-exp(-point_r.^2/2))./(2*pi*A*point_r.^2);
    points(:,1)=Sv*chi./(A*(1+chi.^2));
    points(:,2)=-Sv./(A*(1+chi.^2));
    stability=zeros(size(point_r));
    
    for j=1:numel(point_r)
        fi_war=fi(point_r(j));
        if fi_war<0
            war=A/(St*abs(fi_war)); % 9.01.2019: There was a mistake here before with -1 at the end!
        elseif fi_war>0
            war=A/sqrt(fi_war);
        end
        if war>=1
            stability(j)=1;
        else
            stability(j)=2;
        end
    end
end

function wartosc=krzywa(r,A,Sv)
wartosc=r.*sqrt(1+((1-exp(-r.^2/2))./(2*pi*A*r.^2)).^2)-Sv/A;
end

function wart=fi(r)
wart=(1-exp(-r.^2/2)).*((1-exp(-r.^2/2))/r.^2-exp(-r.^2/2))./(r*(2*pi)).^2;
end

