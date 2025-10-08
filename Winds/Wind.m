classdef Wind < matlab.mixin.Heterogeneous
    %WIND Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Abstract)
        Name string
        wind_model_type (1, 1) {mustBeA(wind_model_type, "WindModelType")}
        include_turbulence (1, 1) logical 
        fields_to_copy (:, 1) string % Fields to copy to the model
    end
    
   methods(Sealed)
        % This section is for methods that call the abstract methods of
        % arrays of subclasses.
        % Example:
        % function constraint_values = evaluate_constraint(aryIn, outputs)
        % 
        %     n = numel(aryIn);
        %     constraint_values = zeros(size(aryIn));
        % 
        %     for c = 1:n
        %         constraint_values(c) = evaluate(aryIn(c), outputs);
        %     end
        % end
    end

    methods(Abstract)
        % Abstract methods that every subclass must implement
        %evaluate(obj, heights) %- implement?

        save(objIn)

        obj = modify(obj, modified_properties)
    end
end

