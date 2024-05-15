function [pasi_err] = absolute_pasi_error(time, fun, time_sim, fun_sim, tot_H, tot_P, PASI_BASELINE)
pasi_err = [];
for i=1:length(time)
    if(isnan(fun(i)))
        continue;
    end
    for j=1:length(time_sim)
        if(time_sim(j)>time(i))
            pasi_err = [pasi_err PASI_BASELINE*(fun(i)-fun_sim(j))/(tot_P-tot_H)];
            break;
        end
    end
end
end