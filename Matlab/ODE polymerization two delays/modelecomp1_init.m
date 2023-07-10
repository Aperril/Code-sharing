function [Xl] = modelecomp1_init(t,X,param)
%Modèle 2 initial pour le complexe
d1=param(1);
d2=param(2);
g1=param(3);
g2=param(4);

Xl=[0;0;0;0;0;0];
%[pc,Ua1,Ua2,psc,C]

Xl(1)=-d1*X(1)*X(2)-d2*X(1)*X(3);
Xl(2)=-d1*X(1)*X(2); %+g1*X(4);
Xl(3)=-d2*X(1)*X(3); %+g2*X(5);
%Xl(4)=-g1*X(4);
%Xl(5)=-g2*X(5);
Xl(4)=0;
Xl(5)=d1*X(1)*X(2);
Xl(6)=d2*X(1)*X(3);
end

