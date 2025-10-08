classdef LaunchSite
    %LAUNCHSITE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name string
        
        latitude double % [deg]
        longitude double % [deg]
        elevation_angle double % [deg] Angle from horizontal
        compass_angle double % [deg] Angle from North (+90 degrees is East)
        elevation double % [m] Initial altitude above MSE
        rail_length double % [m]
        rail_mu double = 0.4 % Coefficient of friction
    end
    
    methods
        function obj = LaunchSite(Name, latitude, longitude, elevation_angle, compass_angle, elevation, rail_length)
            %LAUNCHSITE Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.Name = Name;

            obj.latitude = latitude;
            obj.longitude = longitude;
            obj.elevation_angle = elevation_angle;
            obj.compass_angle = compass_angle;
            obj.elevation = elevation;
            obj.rail_length = rail_length;
        end
        
        function save(obj)
            %SAVE Save launch site to 'LaunchSites' folder in a convenient way.
            
            launch_site = obj;
            
            path = string("LaunchSites\" + obj.Name + ".mat");

            save(path, "launch_site")
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

