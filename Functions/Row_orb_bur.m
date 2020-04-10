function [ r0 ] = Row_orb_bur( m )
%ROWNANIE_ORBITY wyznacza promień orbity stacjonarnej wiru Burgersa bez
%wpływu grawitacji biorąc parametr m=A/St.


syms r
assume(r,{'positive','real'})
r0_sym=solve((m)^(1/2)*r^2-(1-exp(-r^2/2))/(2*pi)==0);
if isempty(r0_sym)==1
    r0=NaN;
else
    r0=double(r0_sym);
end