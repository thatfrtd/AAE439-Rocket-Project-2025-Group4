classdef Motor < matlab.mixin.Heterogeneous
    %MOTOR Summary of this class goes here
    %   Detailed explanation goes here
    % TODO:
    % MAKE ABSTRACT CLASS AND HAVE SUBCLASSES FOR SOLID, LIQUID, etc
    
    properties(Abstract)
        Name string
        
        thrust_table (:, 2) double
        mass_table (:, 2) double
        moi_table (:, 3, 3) double
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
        save(objIn)

        obj = modify(obj, modified_properties)
    end

    methods
        function obj = set_thrust_table_from_file(obj, file, time_column, ...
                thrust_column, convert_to_metric)
            %SET_THRUST_TABLE_FROM_FILE Read in thrust data
            %   File must be compatible with 'readmatrix' PUT IN MOTOR CLASS
            
            Tx_import = readmatrix(file);
            obj.thrust_table = zeros([length(Tx_import), 2]);
            obj.thrust_table(:, 1) = Tx_import(1:length(Tx_import), time_column);
            obj.thrust_table(:, 2) = Tx_import(1:length(Tx_import), thrust_column);

            if convert_to_metric
                obj.thrust_table(:, 2) = obj.Tx(:, 2) .* 4.4482216; % lbf to N
            end
        end


        function obj = set_thrust_mass_tables_from_file(obj, file, time_column, ...
                thrust_column, mass_column, convert_to_metric)
            %SET_THRUST_MASS_TABLES_FROM_FILE Read in thrust and mass data
            %   File must be compatible with 'readmatrix'
            
            Tx_import = readmatrix(file);
            obj.thrust_table = zeros([length(Tx_import), 2]);
            obj.mass_table = zeros([length(Tx_import), 2]);
            obj.thrust_table(:, 1) = Tx_import(1:length(Tx_import), time_column);
            obj.thrust_table(:, 2) = Tx_import(1:length(Tx_import), thrust_column);
            obj.mass_table(:, 1) = Tx_import(1:length(Tx_import), time_column);
            obj.mass_table(:, 2) = Tx_import(1:length(Tx_import), mass_column);

            if convert_to_metric
                obj.thrust_table(:, 2) = obj.thrust_table(:, 2) .* 4.4482216; % lbf to N
                obj.mass_table(:, 2) = obj.mass_table(:, 2) .* 0.454; % lb to kg
            end
        end
    end
end

