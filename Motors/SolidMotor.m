classdef SolidMotor < Motor
    %SolidMotor Summary of this class goes here
    %   Detailed explanation goes here

    properties
        % Inherited Properties
        Name
        
        thrust_table
        mass_table
        moi_table

        % Class Properties
        dry_mass double % [kg]
        total_impulse double
        total_prop_mass double % [kg]

        motor_length double % [m]
        motor_diameter double % [m] 
        prop_radius double % [m]

        bulkhead_m double % [kg]
        casing_m double % [kg]
        liner_m double % [kg]
        nozzle_m double % [kg]

        bulkhead_r double % [m]
        casing_r double % [m]
        liner_r double % [m]
        nozzle_r double % [m]
        prop_r double % [m]
    end

    methods
        function obj = SolidMotor(Name, thrust_table, mass_table, moi_table, dry_mass, total_impulse, total_prop_mass, motor_length, motor_diameter, prop_radius, bulkhead_m, casing_m, liner_m, nozzle_m, bulkhead_r, options)
            %MOTOR Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                Name

                thrust_table
                mass_table = []
                moi_table = []

                dry_mass = []
                total_impulse = []
                total_prop_mass = []
        
                motor_length = []
                motor_diameter = []
                prop_radius = []
        
                bulkhead_m = []
                casing_m = []
                liner_m = []
                nozzle_m = []

                bulkhead_r = []
                options.casing_r = bulkhead_r + motor_length / 2
                options.liner_r = bulkhead_r + motor_length / 2
                options.nozzle_r = motor_length - bulkhead_r
                options.prop_r = bulkhead_r + motor_length / 2
            end
            
            if nargin > 0
                obj.Name = Name;
                obj.thrust_table = thrust_table;
                obj.mass_table = mass_table;
                obj.moi_table = moi_table;
                obj.dry_mass = dry_mass;
                obj.total_impulse = total_impulse;
                obj.total_prop_mass = total_prop_mass;
                obj.motor_length = motor_length;
                obj.motor_diameter = motor_diameter;
                obj.prop_radius = prop_radius;
                obj.bulkhead_m = bulkhead_m;
                obj.casing_m = casing_m;
                obj.liner_m = liner_m;
                obj.nozzle_m = nozzle_m;
                obj.bulkhead_r = bulkhead_r;
                obj.casing_r = options.casing_r;
                obj.liner_r = options.liner_r;
                obj.nozzle_r = options.nozzle_r;
                obj.prop_r = options.prop_r;
            end
        end

        function read_motor_file(obj)
            %READ_MOTOR Read in motor file from thrustcurve.org
        end

        function save(obj)
            %SAVE Save motor to 'Motors' folder in a convenient way.
            
            motor = obj;
            
            path = string("Motors\" + obj.Name + ".mat");

            save(path, "motor")
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

