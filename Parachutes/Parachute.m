classdef Parachute
    %PARACHUTE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name string

        CD double
        radius double % [m]
        deploy_height double = 0 % [m]
    end
    
    methods
        function obj = Parachute(Name, CD, radius, deploy_height)
            %PARACHUTE Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                Name
                CD
                radius
                deploy_height = 0
            end
           
            obj.Name = Name;
            obj.CD = CD;
            obj.radius = radius;
            obj.deploy_height = deploy_height;
        end
        
        function save(obj)
            %SAVE Save parachute to 'Parachutes' folder in a convenient way.
            
            parachute = obj;
            
            path = string("Parachutes\" + obj.Name + ".mat");

            save(path, "parachute")
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

