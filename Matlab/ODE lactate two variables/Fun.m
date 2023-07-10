function[Ft] = F(t,param)
ti=param(7);
tf=param(8);
alpf=param(9);
F0=param(10);

tend=t(end);
Ft=F0*ones(size(t));

for i=1:size(t,2)
        if t(i)<=ti
            Ft(i)=F0;
        elseif (t(i)>ti) && (t(i)<tf)
            Ft(i)=t(i)*F0*alpf/(tf-ti)+(tf*F0-ti*F0*(1+alpf))/(tf-ti);
        else
            Ft(i)=(1+alpf)*F0;
        end
end
end


