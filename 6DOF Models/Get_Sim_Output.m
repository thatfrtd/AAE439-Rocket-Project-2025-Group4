function [outputs] = Get_Sim_Output(simOut, requested_outputs)
%GET_SIM_OUTPUT Extract outputs from a Simulink output object
%   Convenience function for getting a struct of simulation outputs from a
%   Simulink SimulationOutput object. If no outputs are given for
%   requested_outputs, it returns all of the output variables.
arguments
    simOut 
    requested_outputs = string(getElementNames(simOut.yout))
end

outputs = struct();

% Extract simulation outputs
outputs.time = simOut.tout;
for o = 1:length(requested_outputs)
    name = requested_outputs(o);
    
    if name ~= ""
        outputs.(name) = squeeze(simOut.yout.get(name).Values.Data);
    end
end
end

