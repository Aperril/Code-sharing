%%%%%%%%%%% MODELE POUR CROISSANCE GLIOME%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

%%%%%%%%%%%LES PARAMETRES%%%%%%%%%
tspan=[0,200];        %intervalle de temps
Lr=0.313; Gr=5; Or=8.34;  %concentrations arterielles

%consommation cellulaire
%Ta=0; Tn=0; Td=0; T=[Ta;Tn;Td]; %PP
%ba=0; bn=0; bd=0; b=[ba;bn;bd]; %PP
%ga=0; gn=0; gd=0; g=[ga;gn;gd]; %PP

Ta=0.6; Tn=0.5; Td=1; %Ta=0.05; Tn=0.05; Td=0.195; 
T=[Ta;Tn;Td]; %PP
ba=0.001; bn=0.0005; bd=0.001; %0.005; bn=0.005; bd=0.05; 
b=[ba;bn;bd]; %PP
ga=0.15; gn=0.08; gd=0.05; %0.9; gn=0.9; gd=0.8; 
g=[ga;gn;gd]; %PP

%michaelis-mentens
SH=0.0361*8.6; 
c=20.14; u=[SH;c];
pe=3; pn=6; pa=3; pd=3*pa; pc=8; p=[pa;pn;pd;pc;pe]; %PP pd
ke=2.5; ka=2.5; kn=2.5; kd=3*ka; kc=0.95; k=[ka;kn;kd;kc;ke]; %PP kd

%echanges max
co=1;
m1=0.217; m2=0.747; m3=co*m2; m=[m1;m2;m3]; %PP m3
%m1=0.217; m2=0.500; m3=m2; m=[m1;m2;m3]; %PP m3 %Essais
h1=0.041; h2=0.147; h3=co*h2; h4=0.239; h=[h1;h2;h3;h4]; %PP h3
%h1=0.00523; h2=0.00102; h3=h2; h4=0.118; h=[h1;h2;h3;h4]; %PARAM AUBERT PPh3
%a1=0.243; a3=0.50; a2=co*a3; a4=0.1061; a5=co*a4; a6=0.25; a7=0.00243; a8=co*a7; %PP a2 a5 a8 / a2 a3
a1=24.3; a3=0.50; a2=a3; a4=106.1; a5=co*a4; a6=0.25; a7=0.00243; a8=co*a7; %AUB
a=[a1;a2;a3;a4;a5;a6;a7;a8];

%ajust
e1=0.01; e2=0.0001;
%an=0.00001; aa=0.00001; ad=0.01; %PP
an=1; aa=1; ad=0.00001; %PP
v=0.0000001;

%la fonction F
F1=100; F2=10; F3=0.036; F4=0.048; F=[F1;F2;F3;F4]; %Per1, Per2, Plat1, Plat2

%les params
X0=ones(17,1);
X0(1)=0.026; X0(2)=0.026; X0(3)=X0(2); X0(4)=7; 
X0(5)=1.2; X0(6)=1.2; X0(7)=X0(6); X0(8)=1.2; X0(9)=4.56;
X0(10)=1; X0(11)=1.09; X0(12)=X0(11); X0(13)=1; X0(14)=0.35;
X0(15)=3; X0(16)=2*X0(15); X0(17)=0.1;

%Tumeur à 0
%X0(3)=0; X0(7)=0; X0(12)=0; X0(17)=0;
%m3=0; h3=0; a2=0; a5=0; a8=0; m=[m1;m2;m3]; h=[h1;h2;h3;h4]; a=[a1;a2;a3;a4;a5;a6;a7;a8];

param=[Lr;Gr;Or;T;b;g;u;p;k;m;h;a;e1;e2;an;aa;ad;F;v];

%%%%%%ESSAI F%%%%%%%
%figure(1)
%t_test=linspace(tspan(1),tspan(2),10000);
%plot(t_test,Fun(t_test,param),'b')
%xlabel('temps')
%ylabel('impulsion F')
%title('Representation de F')
%set(gca,'fontsize',16)

%%%%%%%DYNAMIQUE%%%%%%
option_ode=[];
[t,Xl]=ode45(@modeleout,tspan,X0,option_ode,param);

%%%%%%AFFICHAGE%%%%%%
figure(2)
subplot(2,2,1)
hold on;
plot(t,Xl(:,1),'b','linewidth',2)
plot(t,Xl(:,2),'k','linewidth',2)
plot(t,Xl(:,3),'r','linewidth',2)
plot(t,Xl(:,4),'g','linewidth',2)
xlabel('time (s)')
ylabel('concentration (mmol/L)')
title('Oxygen')
legend('neuronal','astrocytic','tumor','capillary')
set(gca,'fontsize',16)

subplot(2,2,2)
hold on;
plot(t,Xl(:,5),'b','linewidth',2)
plot(t,Xl(:,6),'k','linewidth',2)
plot(t,Xl(:,7),'r','linewidth',2)
plot(t,Xl(:,8),'m','linewidth',2)
plot(t,Xl(:,9),'g','linewidth',2)
xlabel('time (s)')
ylabel('concentration (mmol/L)')
title('Glucose')
legend('neuronal','astrocytic','tumor','extracell','capillary')
set(gca,'fontsize',16)

subplot(2,2,3)
hold on;
plot(t,Xl(:,10),'b','linewidth',2)
plot(t,Xl(:,11),'k','linewidth',2)
plot(t,Xl(:,12),'r','linewidth',2)
plot(t,Xl(:,13),'m','linewidth',2)
plot(t,Xl(:,14),'g','linewidth',2)
xlabel('time (s)')
ylabel('concentration (mmol/L)')
title('Lactate')
legend('neuronal','astrocytic','tumor','extracell','capillary')
set(gca,'fontsize',16)

Qe=(Xl(:,15)+Xl(:,16)+Xl(:,17))./100;

subplot(2,2,4)
hold on;
plot(t,Xl(:,15)./Qe,'b','linewidth',2)
plot(t,Xl(:,16)./Qe,'k','linewidth',2)
plot(t,Xl(:,17)./Qe,'r','linewidth',2)
xlabel('time (s)')
ylabel('% cells')
title('Cells')
legend('neuronal','astrocytic','tumor')
set(gca,'fontsize',16)

%Necrose ?
figure(3)
subplot(2,2,1)
hold on;
plot(t,Xl(:,15)./Qe,'b','linewidth',2)
plot(t,Xl(:,16)./Qe,'k','linewidth',2)
plot(t,Xl(:,17)./Qe,'r','linewidth',2)
xlabel('temps')
ylabel('concentration')
title('Cellules')
legend('neuronal','astrocytaire','tumoral')
set(gca,'fontsize',16)

N=ones(size(Xl(:,17)));
N(1)=0;
for i=2:size(Xl(:,17),1)
   N(i)=N(i-1)+max(Xl(i-1,17)-Xl(i,17),0);
end

subplot(2,2,2)
hold on;
plot(t,Xl(:,17),'r','linewidth',2)
plot(t,N,'c','linewidth',2)
plot(t,Xl(:,17)+N,'k','linewidth',2)
xlabel('temps')
ylabel('concentration')
title('Cellules')
legend('tumoral','necrose', 'tumeur totale')
set(gca,'fontsize',16)

%PIES
f=size(t,1);
t1=1;
t2=floor(f/3);
t3=floor(2*f/3);
t4=f;

Q1=(Xl(t1,15)+Xl(t1,16)+Xl(t1,17))/100;
X1=[Xl(t1,15)/Q1; Xl(t1,16)/Q1; Xl(t1,17)/Q1];
figure(4)
subplot(2,2,1)
labels = {'Neurons','Astrocytes','Tumor'};
pie(X1,labels)
title('Cells at t=0s')
colormap([0 0 1; 
          .5 .5 .5; 
          1 0 0])
          %0 1 0  

Q2=(Xl(t2,15)+Xl(t2,16)+Xl(t2,17))/100;
X2=[Xl(t2,15)/Q2; Xl(t2,16)/Q2; Xl(t2,17)/Q2];
subplot(2,2,2)
labels = {'Neurons','Astrocytes','Tumor'};
pie(X2,labels)
title('Cells at t=tf/3s')

Q3=(Xl(t3,15)+Xl(t3,16)+Xl(t3,17))/100;
X3=[Xl(t3,15)/Q3; Xl(t3,16)/Q3; Xl(t3,17)/Q3];
subplot(2,2,3)
labels = {'Neurons','Astrocytes','Tumor'};
pie(X3,labels)
title('Cells at t=2tf/3s')

Q4=(Xl(t4,15)+Xl(t4,16)+Xl(t4,17))/100;
X4=[Xl(t4,15)/Q4; Xl(t4,16)/Q4; Xl(t4,17)/Q4];
subplot(2,2,4)
labels = {'Neurons','Astrocytes','Tumor'};
pie(X4,labels)
title('Cells at tf=200000s')