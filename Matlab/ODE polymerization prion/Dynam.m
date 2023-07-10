function [comp,XLs] = Dynam(PARAM,t_ex) %=[lags,m,p,param]
global sol_init lags

%paramétres
param=PARAM(4:end); lags=PARAM(1); m=PARAM(2); p=PARAM(3);
xf=floor(param(1));x0=floor(param(2)); param(1)=xf; param(2)=x0;
xs=floor(param(xf+14));%aj entier
option_ode=[];
T=t_ex(end);

%Initialisation
X0_init=zeros(2*xf+x0+xs+2,1);
X0_init(2*xf+x0+3)=m; %monomères sains
%X0_init(1)=m; %monomères altérés
X0_init(x0+2*xf)=p; %prions sains

%Dynamique
tspan_init=[0,lags];
sol_init=ode45(@modeledis_init,tspan_init,X0_init,option_ode,param);
tspan=[0,T];
sol_fin=dde23('modeledis',lags,'hist',tspan,option_ode,param);

%Estimation sortie
XLs=ones(2*xf+x0+xs+2,size(t_ex,1));
comp=0;

for i=1:size(t_ex,1);
if t_ex(i)<lags
    XLs(:,i)=deval(sol_init,t_ex(i)); comp=comp+1;
else
   XLs(:,i)=deval(sol_fin,t_ex(i));
end
end
end