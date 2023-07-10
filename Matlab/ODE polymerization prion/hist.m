function [y] = hist(t,param)
global sol_init lags
y=deval(sol_init,t+lags);
end