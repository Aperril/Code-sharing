function [Xl] = modeleout(t,X,param)
Lr=param(1); Gr=param(2); Or=param(3);
Ta=param(4); Tn=param(5); Td=param(6);
ba=param(7); bn=param(8); bd=param(9);
ga=param(10); gn=param(11); gd=param(12);
SH=param(13); cu=param(14);
pa=param(15); pn=param(16); pd=param(17); pc=param(18); pe=param(19);
ka=param(20); kn=param(21); kd=param(22); kc=param(23); ke=param(24);
m1=param(25); m2=param(26); m3=param(27);
h1=param(28); h2=param(29); h3=param(30); h4=param(31);
a1=param(32); a2=param(33); a3=param(34); a4=param(35);
a5=param(36); a6=param(37); a7=param(38); a8=param(39);
e1=param(40); e2=param(41);
an=param(42); aa=param(43); ad=param(44);
v=param(49);

Xl=zeros(17,1);

%Oxygene
Xl(1)=m2*(X(4)*SH/(X(4)+cu)-X(1))-6*Tn*X(1)*X(5)-6*gn*X(1)*X(10);
Xl(2)=m1*(X(4)*SH/(X(4)+cu)-X(2))-6*Ta*X(2)*X(6)-6*ga*X(2)*X(11);
Xl(3)=m3*(X(4)*SH/(X(4)+cu)-X(3))-6*Td*X(3)*X(7)-6*gd*X(3)*X(12);
Xl(4)=1/e1*(m2*X(1)+m1*X(2)+m3*X(3)-SH*(m1+m2+m3)*X(4)/(X(4)+cu)+Fun(t,param)*(Or-X(4)));

%Glucose
Xl(5)=h1*(X(8)/(X(8)+pe)-X(5)/(X(5)+pn))-Tn*X(1)*X(5)-bn*X(5);
Xl(6)=h2*(X(8)/(X(8)+pe)-X(6)/(X(6)+pa))-Ta*X(2)*X(6)-ba*X(6);
Xl(7)=h3*(X(8)/(X(8)+pe)-X(7)/(X(7)+pd))-Td*X(3)*X(7)-bd*X(7);
Xl(8)=h1*(X(5)/(X(5)+pn)-X(8)/(X(8)+pe))+h2*(X(6)/(X(6)+pa)-X(8)/(X(8)+pe))+h3*(X(7)/(X(7)+pd)-X(8)/(X(8)+pe))+h4*(X(9)/(X(9)+pc)-X(8)/(X(8)+pe));
Xl(9)=1/e1*(h4*(X(8)/(X(8)+pe)-X(9)/(X(9)+pc))+Fun(t,param)*(Gr-X(9)));

%Lactates
Xl(10)=a1*(X(13)/(X(13)+ke)-X(10)/(X(10)+kn))+a2*(X(12)/(X(12)+kd)-X(10)/(X(10)+kn))+a3*(X(11)/(X(11)+ka)-X(10)/(X(10)+kn))-gn*X(10)*X(1)+2*bn*X(5);
Xl(11)=a4*(X(13)/(X(13)+ke)-X(11)/(X(11)+ka))+a3*(X(10)/(X(10)+kn)-X(11)/(X(11)+ka))+a7*(X(14)/(X(14)+kc)-X(11)/(X(11)+ka))-ga*X(11)*X(2)+2*ba*X(6);
Xl(12)=a2*(X(10)/(X(10)+kn)-X(12)/(X(12)+kd))+a5*(X(13)/(X(13)+ke)-X(12)/(X(12)+kd))+a8*(X(14)/(X(14)+kc)-X(12)/(X(12)+kd))-gd*X(12)*X(3)+2*bd*X(7);
Xl(13)=a1*(X(10)/(X(10)+kn)-X(13)/(X(13)+ke))+a5*(X(12)/(X(12)+kd)-X(13)/(X(13)+ke))+a6*(X(14)/(X(14)+kc)-X(13)/(X(13)+ke))+a4*(X(11)/(X(11)+ka)-X(13)/(X(13)+ke));
Xl(14)=1/e1*(a7*(X(11)/(X(11)+ka)-X(14)/(X(14)+kc))+a8*(X(12)/(X(12)+kd)-X(14)/(X(14)+kc))+a6*(X(13)/(X(13)+ke)-X(14)/(X(14)+kc))+Fun(t,param)*(Lr-X(14)));

%Cellules
Xl(15)=an*(1-X(15)/(e2+1*(4*bn*X(1)+26*gn*X(10)*X(5)+30*Tn*X(1)*X(5))))*X(15)-v*X(15)*X(17); %coeff 1 usage
Xl(16)=aa*(1-X(16)/(e2+1*(4*ba*X(2)+26*ga*X(11)*X(6)+30*Ta*X(2)*X(6))))*X(16)-v*X(16)*X(17); %coeff 1 usage
Xl(17)=ad*(1-X(17)/(e2+4*bd*X(3)+26*gd*X(12)*X(7)+30*Td*X(3)*X(7)))*X(17);
end

