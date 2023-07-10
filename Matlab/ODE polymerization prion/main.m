%%%%%%%%%%Essai de dynamiques
clear all, close all, clc,
warning ('off','all');
disp(['Voici la dynamique globale pour un peptide Abeta pouvant former']);
disp(['un complexe avec le prion sain.']);
disp([' ']);

%Paramètres
T=linspace(0,120,5000)';                   % Intervalle temps évaluation
xf=40; x0=20; xs=13;                       % x critique
b=0.01; ba=0.1; bf=0.1;                    % taux dépolym
del=0.8; gam=0.2; gamf=1;                % complexe et nettoyage
eps=0.002; alp=0.01;                       % faible phénomène
am=0.2; bm=0.9; ae=0.9; be=0.8;              % micellisation
a=0.08*ones(xf-1,1); %0.1                   % taux polym
%a(1:10)=0.09;                             % taux polym altéré
param=[xf;x0;b;ba;bf;del;gam;gamf;alp;eps;a;am;bm;ae;be;xs];

%Dynamique
lags=1;                                  % retard
m=50;                                     % monomères sains initiaux
p=50;                                     % prions sains initiaux
PARAM=[lags;m;p;param];
[comp,XL]=Dynam(PARAM,T);

%Affichage
figure(1)
subplot(2,2,3)
hold on;
plot(T,XL(1,:),'b','linewidth',2)
plot(T,x0*XL(x0,:),'r','linewidth',2)
plot(T,x0*XL(x0+2*xf+2,:),'g','linewidth',2)
plot(T,x0*XL(x0+1,:),'m','linewidth',2)
%plot(T,Q,'y')
xlabel('temps (t)')
ylabel('masse (g)')
title('Trajectoires : polymères')
legend('monomères','oligomères','complexes','o-plaques') %,'masse Q')
xm=max([max(XL(1,:)),x0*max(XL(x0,:)),x0*max(XL(x0+2*xf+2,:)),x0*max(XL(x0+1,:))]);
line([lags lags], [0 xm],'Color', 'k','Linestyle',':');
set(gca,'fontsize',12)
axis([0 T(end) 0 xm])
disp(['Figure 1.1 : monomères et oligomères']);

subplot(2,2,2)
hold on;
plot(T,(x0-2)*XL((x0-2),:),'k','linewidth',2)
plot(T,2*XL((2),:),'m','linewidth',2)
plot(T,x0*XL(x0,:),'b','linewidth',2)
%plot(T,Q,'y')
xlabel('temps (t)')
ylabel('masse (g)')
title('Trajectoires : polymères proto')
legend('proto oligo x0-2','proto-oligo 2','oligomères x0') %,'masse Q')
xm=max([x0*max(XL(x0,:)),(x0-2)*max(XL(x0-2,:)),(2)*max(XL(2,:))]);
line([lags lags], [0 xm],'Color', 'k','Linestyle',':');
set(gca,'fontsize',12)
axis([0 T(end) 0 xm])
disp(['Figure 1.2 : quelques étapes proto-oligomériques']);

subplot(2,2,1)
hold on;
plot(T,XL(x0+2*xf,:),'b','linewidth',2)
plot(T,XL(x0+2*xf+1,:),'r','linewidth',2)
plot(T,XL(x0+2*xf+2,:),'g','linewidth',2)
%plot(T,Qp,'y')
xlabel('temps (t)')
ylabel('masse (g)')
title('Trajectoires : prion')
legend('prions sains','prions malades','complexes') %,'masse Qp')
xm=max([max(XL(x0+2*xf,:)),max(XL(x0+2*xf+1,:)),max(XL(x0+2*xf+2,:))]);
line([lags lags], [0 xm],'Color', 'k','Linestyle',':');
set(gca,'fontsize',12)
axis([0 T(end) 0 xm])
disp(['Figure 1.3 : étapes pour le prion']);

subplot(2,2,4)
hold on;
plot(T,2*XL(x0+2,:),'k','linewidth',2)
plot(T,3*XL(x0+3,:),'b','linewidth',2)
plot(T,4*XL(x0+4,:),'r','linewidth',2)
%plot(T,Q,'y')
xlabel('temps (t)')
ylabel('masse (g)')
title('Trajectoires : fibrilles')
legend('fibrilles 2','fibrilles 3','fibrilles 4') %,'masse Q')
xm=max([max(2*XL(x0+2,:)),3*max(XL(x0+3,:)),max(4*XL(x0+4,:))]);
line([lags lags], [0 xm],'Color', 'k','Linestyle',':');
set(gca,'fontsize',12)
axis([0 T(end) 0 xm])
disp(['Figure 1.4 : quelques étapes fibrillaires']);

%Suivi masse renormalisée
Qi=[1:xs]*XL(((x0+2*xf+3):(x0+2*xf+xs+2)),:);
Q=XL(1,:)+[2:x0]*XL((2:x0),:)+x0*XL((x0+1),:)+[2:xf]*XL(((x0+2):(x0+xf)),:)+[2:xf]*XL(((x0+xf+1):(x0+2*xf-1)),:)+x0*XL(x0+2*xf+2,:)+Qi;
Qp=XL(x0+2*xf,:)+XL(x0+2*xf+1,:)+XL(x0+2*xf+2,:);
test=0;
for i=1:size(Q,2)
if round(Q(i)*1000)/1000~=round(Q(1)*1000)/1000; test=1;
end
end
if test==0 disp(['La masse de polymères est restée constante à ' num2str(Q(1)) '.']);
    else disp(['Attention la masse des polymères a été fluctuante !']);
end
disp([' ']);

%1er moment 
figure(2)
subplot(2,2,1)
plot(T,Q,'b','linewidth',2)
xlabel('Temps (t)')
ylabel('1er moment')
title('Evolution du 1er moment')
disp(['Figure 2.1 : évolution du premier moment']);

%1er moment sans <x0
subplot(2,2,2)
Q0=x0*XL(x0,:)+x0*XL((x0+1),:)+[x0:xf]*XL((2*x0:(x0+xf)),:)+[x0:xf]*XL(((2*x0+xf-1):(x0+2*xf-1)),:)+x0*XL(x0+2*xf+2,:);
plot(T,Q0,'b','linewidth',2)
xlabel('Temps (t)')
ylabel('1er moment')
title('Evolution du 1er moment sans <x0')
disp(['Figure 2.2 : évolution du premier moment sans <x0']);

%Coeff thioflav
subplot(2,2,3)
c0=22+(9/25).*Q0;
plot(T,c0,'b','linewidth',2)
xlabel('Temps (t)')
ylabel('1er moment')
title('Evolution du 1er moment avec coeff thioflav')
disp(['Figure 2.3 : 1er moment avec coeff thioflav']);
disp([' ']);