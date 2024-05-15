function [sbml] = sbml_set_parameter_value(sbml, par, value)
%SBML_SET_PARAMETER_VALUE Summary of this function goes here
%   Detailed explanation goes here
    % parameter already exists
    for i = 1:length(sbml.Parameters)
        if(strcmp(sbml.Parameters(i).Name, par))
            sbml.Parameters(i).Value = value;
            return;
        end
    end
    % parameter does not exist
    sbml.addparameter(par, value);	
end

