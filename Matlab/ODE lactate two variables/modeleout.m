function [Xl] = modeleout(t,X,param)
%modele general pour les lactates

kap=param(1); k=param(2); kp=param(3); eps=param(4);
L=param(5); Cj=param(6);

Xl=[0;0];

J=Cj/(X(1)+eps);

Xl(1)=J+kap*(X(2)/(kp+X(2))-X(1)/(k+X(1)));
Xl(2)=(1/eps)*(Fun(t,param)*(L-X(2))+kap*(X(1)/(k+X(1))-X(2)/(kp+X(2))));
end

