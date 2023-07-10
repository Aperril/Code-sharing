function [y] = sig(params,t)
a=params(1); M=params(2); P=params(3); T0=params(4); k=params(5);

for i=1:size(t,1)
dem=1+(M/P-1)*exp(-k*(t(i)-T0));
y(i)=a+M/dem;
end

y=y';
end

