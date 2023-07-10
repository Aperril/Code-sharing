function [y] = sig(params,t)
A1=params(1); A2=params(2); k=params(3); t0=params(4);

for i=1:size(t,1)
y(i)=(A2-A1)/(1+exp(k*(t(i)-t0)))+A2;
end
y=y';
end

