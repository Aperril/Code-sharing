%%%%%Modèle avec deux retards (1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc
global sol_init lags2 lags1 sol_med 

%Choisis tes paramètres de dynamique
d1=0.2; g1=0.4; lags1=2;
d2=0.4; g2=0.6; lags2=1;
X0_init=[20;10;10;0;0;0]; %[pc,Ua1,Ua2,psc,C1,C2]

%formalités
option_ode=[];
lags=[lags1;lags2];
param=[d1;d2;g1;g2];
n=floor(lags1/lags2);
Tf=20;

%Dynamique initiale
%X0_init=[40;10;10;0;0;0;0;0];
tspan_init=[0,lags(2)];
sol_init=ode45(@modelecomp1_init,tspan_init,X0_init,option_ode,param);
t_init=sol_init.x;
Xl_init=sol_init.y'; 

%Dynamique mediane (1 retard)
tspan_med=[lags2,max(2,n)*lags2];
sol_med=dde23('modelecomp1_med',lags2,'hist1',tspan_med,option_ode,param);
t_med=sol_med.x;
Xl_med=sol_med.y;

%Dynamique finale (2 retards)
tspan_fin=[lags(1),Tf];
sol_fin=dde23('modelecomp1_fin',lags,'hist2',tspan_fin,option_ode,param);
t_fin=sol_fin.x;
Xl_fin=sol_fin.y;

%Troncature
t_med2=linspace(lags2,lags1,50);
Xl_med2=hist2(t_med2,param);

%Affichage total
nm=max(X0_init);
T=[t_init,t_med2,t_fin];
XL=[Xl_init',Xl_med2,Xl_fin];

subplot(2,2,1)
hold on;
plot(T,XL(2,:),'b','linewidth',2)
%plot(T,XL(4,:),'r','linewidth',2)
plot(T,XL(5,:),'k','linewidth',2)
xlabel('time')
ylabel('concentration')
title('Oligomer''s 1 evolution (2 varieties)')
legend('Active oligomer','Complex')
line([lags2 lags2], [0 12],'Color', 'k','Linestyle',':');
line([lags1 lags1], [0 12],'Color', 'k','Linestyle',':');
set(gca,'fontsize',12)
axis([0 20 0 10])

subplot(2,2,2)
hold on;
plot(T,XL(3,:),'b','linewidth',2)
%plot(T,XL(5,:),'r','linewidth',2)
plot(T,XL(6,:),'m','linewidth',2)
xlabel('time')
ylabel('concentration')
title('Oligomer''s 2 evolution (2 varieties)')
legend('Active oligomer','Complex')
line([lags2 lags2], [0 12],'Color', 'k','Linestyle',':');
line([lags1 lags1], [0 12],'Color', 'k','Linestyle',':');
set(gca,'fontsize',12)
axis([0 20 0 10])


subplot(2,2,3)
hold on;
plot(T,XL(1,:),'b','linewidth',2)
plot(T,XL(4,:),'r','linewidth',2)
plot(T,XL(5,:),'k','linewidth',2)
plot(T,XL(6,:),'m','linewidth',2)
xlabel('time')
ylabel('concentration')
title('Prion''s evolution (2 varieties)')
legend('Healthy prion','Altered prion','Complex 1','Complex 2')
line([lags2 lags2], [0 nm],'Color', 'k','Linestyle',':');
line([lags1 lags1], [0 nm],'Color', 'k','Linestyle',':');
set(gca,'fontsize',12)
axis([0 20 0 20.001])