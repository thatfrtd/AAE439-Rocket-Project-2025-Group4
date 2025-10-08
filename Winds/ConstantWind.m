classdef ConstantWind < Wind
    %CONSTANTWIND Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Inherited Properties
        Name
        wind_model_type = WindModelType.Wind_Table
        include_turbulence = false
        fields_to_copy = ["wind_model_type", "include_turbulence", "wind_direction_table", "wind_speed_table"]

        % Class Properties
        wind_direction double
        wind_speed double
    end

    properties (Dependent)
        wind_direction_table % 0 is North, +90 is East
        wind_speed_table
    end

    methods
        function obj = ConstantWind(Name, wind_direction, wind_speed)
            %CONSTANTWIND Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.Name = Name;

            obj.wind_direction = wind_direction;
            obj.wind_speed = wind_speed;
        end

        function wind_direction_table = get.wind_direction_table(obj)            
            time = [0; 1];
            wind_direction_table = [time, [obj.wind_direction; obj.wind_direction]];
        end

        function wind_speed_table = get.wind_speed_table(obj)            
            time = [0; 1];
            wind_speed_table = [time, [obj.wind_speed; obj.wind_speed]];
        end

        function save(obj)
            %SAVE Save wind to 'Winds' folder in a convenient way.
            
            wind = obj;
            
            path = string("Winds\" + obj.Name + ".mat");

            save(path, "wind")
        end

        function obj = modify(obj, modified_properties)
            %MODIFY Modify properties of this class. Must be a table.

            for p = 1:size(modified_properties, 2)
                name = string(modified_properties.Properties.VariableNames(p));
                
                % Perform checks and make corrections if necessary
                
                % Modify
                obj.(name) = modified_properties.(name);
            end
        end
    end
end

