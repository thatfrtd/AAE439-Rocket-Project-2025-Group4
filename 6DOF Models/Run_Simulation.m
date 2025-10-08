function [outputs] = Run_Simulation(model_name, model_inputs, options)
%RUN_SIMULATION Run simulation and return outputs
%   Runs the simulation with inputs in model_inputs object and return requested
%   outputs. Requested outputs must be from output blocks (appear in yout).
%   If the model_inputs object is an array of model_inputs, then
%   simulations for each model_inputs object are run in parallel (if
%   available) and the outputs are returned as an array of outputs.
arguments
    model_name
    model_inputs
    options.stall_time = inf
    options.requested_outputs = []
    options.show_simulation_manager string {mustBeMember(options.show_simulation_manager, ["on", "off"])} = "off"
    options.show_progress string {mustBeMember(options.show_progress, ["on", "off"])} = "off"
end

%% Initialize 6DoF
sim_number = numel(model_inputs);

% Open model so input can be created for it and it can be ran
open(string("6DOF Models/" + model_name));

% Create a simulation input object to store all the input for the simulation
simIn(1:sim_number) = Simulink.SimulationInput(model_name);

% Send data in model_input struct to simulation
input_field_names = fieldnames(model_inputs(1));
for s = 1:sim_number
    for i = 1:length(input_field_names)
        name = string(input_field_names(i));
    
        simIn(s) = setVariable(simIn(s), name, model_inputs(s).(name), Workspace = model_name);
    end
end

%% Run 6DoF
simIn = setModelParameter(simIn, "TimeOut", options.stall_time);

simIn = setModelParameter(simIn, 'SimulationMode', 'accelerator');

% Run sim
w = warning('off','all');
if sim_number > 1
    simOut = parsim(simIn, 'ShowSimulationManager', options.show_simulation_manager, ...
        'ShowProgress', options.show_progress, 'StopOnError', 'off');
else
    simOut = sim(simIn, 'ShowSimulationManager', options.show_simulation_manager, ...
        'ShowProgress', options.show_progress, 'StopOnError', 'off');
end
warning(w);

%% Extract Outputs
if isempty(options.requested_outputs)
    options.requested_outputs = string(getElementNames(simOut(1).yout));
end

outputs(1:sim_number) = struct();
for s = 1:sim_number
    % Extract simulation outputs
    outputs(s).time = simOut(s).tout;
    for o = 1:length(options.requested_outputs)
        name = options.requested_outputs(o);
        
        outputs(s).(name) = squeeze(simOut(s).yout.get(name).Values.Data);
    end
    % Propogate errors
    % ModelStop is the value if stops because Stop block 
    % - Use block path from StopEventSource to find name of stop block? 
    %   Path is just to model...
    outputs(s).Error_Message = simOut(s).ErrorMessage;
    if simOut(s).SimulationMetadata.ExecutionInfo.StopEvent == "Timeout" 
        outputs(s).Timedout = true;
    else
        outputs(s).Timedout = false;
    end
end

end