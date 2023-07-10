function [Xl] = modelelim(t,X,param)
%modele limite eps=0, /!\ seul u retourne
kap=param(1); k=param(2); kp=param(3); eps=param(4); L=param(5);
Cj=param(6);

Xl=[0];

J=Cj/(X(1)+eps);

z=kap*X(1)/(Fun(t,param)*(k+X(1)));
v=(1/2)*(-kp+L-(kap/Fun(t,param))+z+sqrt((kp-L+(kap/Fun(t,param))-z)^2+4*kp*(L+z)));

Xl(1)=J+kap*(v/(kp+v)-X(1)/(k+X(1)));
end

