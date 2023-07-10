function [Xl] = modeledis_init(t,X,param)
%Les parametres
xf=param(1); x0=param(2); xs=param(xf+14);
b=param(3); ba=param(4); bf=param(5);
del=param(6); gam=param(7); gamf=param(8); alp=param(9);
eps=param(10); a=param(11:(xf+9))'; am=param(xf+10); bm=param(xf+11);
ae=param(xf+12); be=param(xf+13);
K=1/(x0-3);

%Les monomères à la polym
rajx0=[2:(x0-1)].*a(2:x0-1);
rajxf=[2:(xf-1)].*a(2:(xf-1));
M=sum(b*X(2:(x0-1)))+sum(bf*X((x0+2):(x0+xf)))+sum(ba*X((x0+xf+1):(x0+2*xf-1)))+b*X(2)+bf*X(x0+2)+ba*X(x0+xf+1);
P=rajx0*X(2:(x0-1))+rajxf*X((x0+2):(x0+xf-1));

%La source à la mic
Mm=sum(bm*X( (x0+2*xf+4) : (x0+2*xf+xs+2) )) +bm*X(x0+2*xf+4);
Pm=sum(am*X( (x0+2*xf+4) : (x0+2*xf+xs+1) ));

%Les monomères a la mic
Me=sum(be*X( (x0+2*xf+4): (x0+2*xf+xs+2) )) +be*X(x0+2*xf+4);
Pe=sum(ae*X( (x0+2*xf+4): (x0+2*xf+xs+1) ));

%La dynamique
Xl=zeros(2*xf+x0+xs+2,1); %[c1,ci,c_{x0-1},c0,ua,cfj,cf0,cai,c0a,pc,psc,Psi,sj]

%X_1=c1
%Xl(1)=-4*1*a(1)*X(1)^2-P*X(1)+M; sans mic
Xl(1)=-4*1*a(1)*X(1)^2-P*X(1)+M-2*ae*X(1)^2-Pe*X(1)+Me;

%X_2 --- X_{x0-2}=ci
for i=2:(x0-2)
    Xl(i)=(i-1)*a(i-1)*X(i-1)*X(1)+b*X(i+1)-X(i)*(b+i*a(i)*X(1))+X(x0)*eps*K;
end

%X_{x0-1} =cx0-1
Xl(x0-1)=(x0-2)*a(x0-2)*X(x0-2)*X(1)-X(x0-1)*((x0-1)*a(x0-1)*X(1)+b); %+X(x0)*eps*K;

%X_{x0}=c0
Xl(x0)=-gam*X(x0)-del*X(x0)*X(x0+2*xf)+(x0-1)*a(x0-1)*X(1)*X(x0-1)-(eps/2)*X(x0);

%X_{x0+1}=ua
Xl(x0+1)=gam*X(x0); 

%X_{x0+2}=cf2
Xl(x0+2)=a(1)*X(1)^2+bf*X(x0+3)-X(x0+2)*(bf+2*a(2)*X(1))-gamf*X(x0+2);

%X_{x0+3} --- X_{x0+xf-1}=cfj
for j=(x0+3):(x0+xf-1);
    Xl(j)=(j-1-x0)*a(j-1-x0)*X(j-1)*X(1)+bf*X(j+1)-X(j)*(bf+(j-x0)*a(j-x0)*X(1))-gamf*X(j);
end

%X_{x0+xf}=cf0
Xl(x0+xf)=-(bf+gamf)*X(x0+xf)+(xf-1)*a(xf-1)*X(1)*X(x0+xf-1);

%X_{x0+xf+1}___X_{x0+2xf-2}=caj
for k=(x0+xf+1):(x0+2*xf-2);
    Xl(k)=ba*X(k+1)-ba*X(k)+gamf*X(k-xf+1);
end

%X_{x0+2xf-1}=ca0
Xl(x0+2*xf-1)=-ba*X(x0+2*xf-1)+gamf*X(x0+xf); 

%X_{x0+2xf}=pc
Xl(x0+2*xf)=-alp*X(x0+2*xf)-del*X(x0+2*xf)*X(x0); 

%X_{x0+2xf+1}=psc
Xl(x0+2*xf+1)=alp*X(x0+2*xf); 

%X_{x0+2xf+2}=Psi
Xl(x0+2*xf+2)=del*X(x0+2*xf)*X(x0);

%miscellisation

%X_{x0+2xf+3}=s1
Xl(x0+2*xf+3)=-2*am*X(x0+2*xf+3)^2-Pm*X(x0+2*xf+3)+Mm;

%X_{x0+2xf+4}=s2
Xl(x0+2*xf+4)=am*X(x0+2*xf+3)^2+bm*X(x0+2*xf+5)-X(x0+2*xf+4)*(am*X(x0+2*xf+3)+bm)+ae*X(1)^2+be*X(x0+2*xf+5)-X(x0+2*xf+4)*(ae*X(1)+be);

%X_{x0+2xf+5}___X_{x0+2xf+xs+1}=sj
for j=(x0+2*xf+5):(x0+2*xf+xs+1);
    Xl(j)=am*X(j-1)*X(x0+2*xf+3)+bm*X(j+1)-X(j)*(am*X(x0+2*xf+3)+bm)+ae*X(j-1)*X(1)+be*X(j+1)-X(j)*(ae*X(1)+be);
end

%X_{x0+2xf+xs+2}=smax
Xl(x0+2*xf+xs+2)=am*X(x0+2*xf+xs+1)*X(x0+2*xf+3)-X(x0+2*xf+xs+2)*bm+ae*X(x0+2*xf+xs+1)*X(1)-X(x0+2*xf+xs+2)*be;
end

