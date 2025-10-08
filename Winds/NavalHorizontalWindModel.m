classdef NavalHorizontalWindModel < Wind
    %NAVALHORIZONTALWINDMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Inherited Properties
        Name
        wind_model_type = WindModelType.HW14_Wind_Model
        include_turbulence = false
        fields_to_copy = ["wind_model_type", "include_turbulence", "day_seconds"]

        % Class Properties
        start_time (1, 1) datetime
    end

    properties (Dependent)
        day_seconds
    end

    methods
        function obj = NavalHorizontalWindModel(Name, start_time, include_turbulence)
            %NAVALHORIZONTALWINDMODEL Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.Name = Name;
            
            obj.include_turbulence = include_turbulence;

            obj.start_time = start_time;
        end

        function day_seconds = get.day_seconds(obj)
            day_of_year =  day(obj.start_time, 'dayofyear');
            second_of_day = obj.start_time.Hour * 3600 + obj.start_time.Minute * 60 + obj.start_time.Second;
            day_seconds = [day_of_year, second_of_day];
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


