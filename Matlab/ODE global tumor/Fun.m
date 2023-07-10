function[Ft] = Fun(t,param)
P1=param(45);
P2=param(46);
F1=param(47);
F2=param(48);

tend=t(end);
Ft=F1*ones(size(t));
n=floor(tend/(P1+P2))+1;

for i=1:size(t,2)
    for j=1:n
        if (t(i)<j*(P1+P2)) && (t(i)>(j-1)*(P1+P2)+P1)
            Ft(i)=F2;
        end
        
    end
end

end


