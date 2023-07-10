function [y] = hist(t,param)
global sol_init sol_med lags2

if t<lags2
    y=deval(sol_init,t);
else
    y=deval(sol_med,t);
end

end