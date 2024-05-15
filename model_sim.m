function [res] = model_sim(sbml,time_bound)
%SBML_SIM_AND_PLOT Summary of this function goes here
%   Detailed explanation goes here
    sbml.getconfigset.StopTime = time_bound;
    res = sbiosimulate(sbml);
end

