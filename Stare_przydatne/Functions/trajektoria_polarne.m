function dx=trajektoria_polarne(t,x,M1r,M2r,M3r)
dx=[x(3);x(4); -M1r*x(1)/2-M2r*x(3)-M3r*sin(x(2))+x(1)*(x(4))^2; M2r*(1-exp(-x(1)^2/2))/(2*pi*x(1)^2)-M3r*cos(x(2))/x(1)-2*x(3)*x(4)/x(1)-M2r*x(4)];
% w tym wzorze jest bład: zamiast + jest * przed ostatnim członen
% pierwszego rownania
end