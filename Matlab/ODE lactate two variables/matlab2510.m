%%%%%%%%%%% MODELE POUR LACTATES DEUX VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

%Présentation
disp(['CINETIQUE DU MODELE SIMPLISTE DE DEUX VARIABLES'])
disp([' '])
disp(['Utilisation d''echanges de type symport.'])
disp(['La fonction F est un unique creneau, J est decroissante en fonction de u'])
disp([' '])

%Mise en place de la dynamique
tspan=[0,10000];        %intervalle de temps
X0=[1.15;1];            %valeurs initiales [u,v]
kap=0.01;               %taux transport max lac cel->cap (MCT) 
k=3.5;                  %cst MM pour u/[H+]
kp=3.5;                 %cst MM pour v/[H+]
eps=0.001;              %mise a l'echelle des eqs
L=0.3;                  %concentraiton lac veineux/arterielle
Cj=5.7*10^(-5);         %pour J
ti=50; tf=100;          %t critiques F
alpf=0.7;               %coeff F
F0=0.012;               %valeurs de base pour F

param=[kap;k;kp;eps;L;Cj;ti;tf;alpf;F0];

figure(1)
subplot(2,1,1)
t0=linspace(0,10000, 600);
plot(t0,Fun(t0,param))
set(gca,'fontsize',16)
xlabel('temps')
ylabel('Blood flow')
title('Fonction F')

subplot(2,1,2)
x0=linspace(0,3,3);
plot(x0,Cj./(x0+eps))
set(gca,'fontsize',16)
xlabel('lactate cellule')
ylabel('Cell flux')
title('Fonction J')
disp(['Figure 1 : F et J'])

%%%%%%%%%%%%%%%%%%%%  MODELE GENERAL  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Dynamique
option_ode=[];
[t,Xl]=ode45(@modeleout,tspan,X0,option_ode,param);

%%%%%%%%%%%%%%%%%%%%  MODELE LIMITE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Dynamique
option_ode=[];
[tl,ul]=ode45(@modelelim,tspan,X0(1),option_ode,param);

vl=ones(size(tl));
for i=1:size(tl,1)
    z=kap*ul(i)/(Fun(tl(i),param)*(k+ul(i)));
    vl(i)=(1/2)*(-kp+L-(kap/Fun(tl(i),param))+z+sqrt((kp-L+(kap/Fun(tl(i),param))-z)^2+4*kp*(L+z)));
end

%%%%%%%%%%%%%%%%%%%%%%%%BORNES%%%%%%%%%%%%%%
%Borne v
Bv0=max(X0(2),(kap+L*(1+alpf)*F0)/F0);
Bv=Bv0*ones(size(Xl,1),1);

%Borne u
J1=Cj*eps;
if J1*(kp+Bv)< kap*kp 
    disp(['La condition de borne est satisfaite.'])
    Bu0=max(X0(1), k*(Bv0/(kp+Bv0)+(J1/kap))/(1-Bv0/(kp+Bv0)-(J1/kap)));
    Bu=Bu0*ones(size(Xl,1),1);
    Bul=Bu0*ones(size(ul,1),1);
else
    disp(['ATTENTION la condition de borne N''est PAS satisfaite !'])
    Bu=0*ones(size(Xl,1),1);
    Bu=0*ones(size(ul,1),1);
end
disp([' '])

%Borne v_lim
Bvl=((kap+L*(1+alpf)*F0)/F0)*ones(1,size(ul,1));

%%%%%%%%%AFFICHAGE%%%%%%%%%%%%%%%%%%%%

%Affichage u,v
disp(['FIGURE 2'])
figure(2)
subplot(2,2,1)
hold on;
plot(t,Xl(:,1),'b')
plot(t,Xl(:,2),'r')
plot(t,Bu,'.')
plot(t,Bv,'m')
xlabel('temps')
ylabel('concentration')
title('Trajectoires des lactates')
legend('lac intracel','lac capillaire', 'borne Li (sous cdt)', 'borne Lc')
set(gca,'fontsize',16)
disp(['Affichage (2,2,1) : Trajectoires lactaires.'])

subplot(2,2,2)
plot(Xl(:,1),Xl(:,2),'b')
xlabel('u')
ylabel('v')
title('Portrait de phase')
set(gca,'fontsize',16)
disp(['Affichage (2,2,2) : Portrait de phase.'])
disp([' '])

%Affichage u',v'
subplot(2,2,3)
hold on;
plot(tl,ul,'b')
plot(tl,vl,'r')
plot(tl,Bul,'.')
plot(tl,Bvl,'m')
xlabel('temps')
ylabel('concentration')
title('Trajectoires des lactates (eps=0)')
legend('lac intracel','lac capillaire' ,'borne Li (sous cdt)', 'borne Lc')
set(gca,'fontsize',16)
disp(['Affichage (2,2,3) : Trajectoires du modèle limite.'])
    
subplot(2,2,4)
plot(ul,vl,'b')
xlabel('u')
ylabel('v')
title('Portrait de phase (eps=0)')
set(gca,'fontsize',16)
disp(['Affichage (2,2,4) : Portrait de phase du modèle limite.'])
disp([' '])