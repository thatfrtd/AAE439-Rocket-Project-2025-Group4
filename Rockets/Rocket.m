classdef Rocket
    %ROCKET Rocket simulation data for 6DoF and DATCOM
    %   Class containing all the parameters that describe a rocket for
    % simulation and putting into Missle DATCOM. Also describes what all
    % the different inputs mean and what their options are.
    %
    % TODO: 
    % Protuberances
    % angled fins with PHIF
    % function to create NACA card
    
    properties
        % General
        Name string = "Rocket"
        mass_dry double % [kg] No motor
        XCG_dry(1, 1) double % [m] Center of gravity measured from nosecone tip
        MOI_dry (3, 3) double % [m2 kg] Moment of inertia with no motor
        LREF(1, 1) double % [m] Reference length
        SREF(1, 1) double % [m2] Reference area
        ROUGH(1, 1) double = 6.35e-6 % [m] Surface Roughness
        DEXIT(1, 1) double = 0 % [m] Exit Diameter, 0 makes base drag from body included in axial force calculations

        % Nosecone
        NC_type string {mustBeMember(NC_type, ["CONICAL", "CONE", "OGIVE", "POWER", "HAACK", "KARMAN"])} = "OGIVE"
            % POWER option needs NC_power
        NC_power(1, 1) double {mustBeInRange(NC_power, 0, 1)} = 0.5 % [] 
            % Exponent, n, for power law shape, (r/R)=(x/L)^n, used for POWER NC_type
            % n = 1: Cone, n = 0.5: Parabola, n=0: Cylinder
        NC_length(1, 1) double % [m] Length                                             
        NC_diam(1, 1) double % [m] Diameter of nosecone  

        % Airframe
        AF1_length(1, 1) double % [m] Airframe length
        AF1_diam(1, 1) double % [m] Airframe diameter

        % SHOULD INCLUDE PROTUBERANCES - launch mount, cameras, etc.
        
        % Fins
            % Array for multiple finsets - need to all have same length
        FN1_number(1, :) {mustBeInteger(FN1_number), mustBeInRange(FN1_number, 1, 8)} = 4 % Number of fins
        FN1_type string {mustBeMember(FN1_type, ["HEX", "ARC", "NACA"])} = "HEX"
            % NACA option needs FN1_NACA
        FN1_NACA string = "" % Control card defining NACA airfoil used for NACA FN1_type 
            % "NACA-{Fin set number (1,2,3,4)}-{Airfoil Type (1,4,5,6,S)}-{Input designation}" 
            % EX: NACA-3-S-1-45.5-6.8
            %       Could create a NACA function that formats stuff for you
            %       but is NACA really useful?
        FN1_root(1, :) double % [m] Root chord length                                     
        FN1_tip(1, :) double % [m] Tip chord length                                       
        FN1_height(1, :) double % [m] Fin span                                            
        FN1_sweep(1, :) double {mustBeInRange(FN1_sweep, 0, 90)} % [degrees] Fin sweep
        FN1_edgeRad(1, :) double = 0.000127 % [m] Fin leading edge radius 
        FN1_XLE(1, :) double % [m] Distance from NC tip to front of leading edge
        
        % Fin Cross-section
            % Array for multiple finsets - need to all have same length
        THICKNESS(1, :) double % [m] Fin thickness ONLY FOR SYMMETRIC FINS (Not for DATCOM - used to calculate ZUPPER and ZLOWER)
            % Following used for HEX and ARC FN1_type
        ZUPPER(1, :) double % [CHORD FRACTION] Thickness of upper surface
        ZLOWER(1, :) double % [CHORD FRACTION] Thickness of lower surface
            % Following used for HEX FN1_type
        LMAXU(1, :) double {mustBeInRange(LMAXU, 0, 0.5)} % [CHORD FRACTION] Leading edge to maximum thickness of upper surface
        LMAXL(1, :) double {mustBeInRange(LMAXL, 0, 0.5)} % [CHORD FRACTION] Leading edge to maximum thickness of lower surface
        LFLATU(1, :) double {mustBeInRange(LFLATU, 0, 1)} % [CHORD FRACTION] Constant thickness section of upper surface
        LFLATL(1, :) double {mustBeInRange(LFLATL, 0, 1)} % [CHORD FRACTION] Constant thickness section of lower surface

        % Tables (Not for DATCOM)
        m_table(:, 2) double % [kg] Mass (column 1 is time) - Same times as thrust - MOVE INTO DEPENDENT?
        COM_table(:, 2) double % [m] Center of mass (column 1 is time) - Same times as thrust
        MOI_table(:, 3, 3) double % [m2 kg] Moment of Inertia (same times as COM_table)
        
        motor {mustBeA(motor, "Motor")} = Motor.empty
        main_parachute {mustBeA(main_parachute, "Parachute")} = Parachute("NoParachute", 0, 0)
        drogue_parachute {mustBeA(drogue_parachute, "Parachute")} = Parachute("NoParachute", 0, 0)
    end

    properties (Dependent)
        Tx % [N] Thrust (column 1 is time)
    end
    
    methods
        function obj = Rocket(mass_dry, XCG_dry, AF1_length, AF1_diam, FN1_root, FN1_tip, FN1_height, FN1_sweep, THICKNESS, LMAXU, NC_length, NC_diam, ZUPPER, options)
            %ROCKET Construct an instance of this class
            %   If only required inputs are given, assumes constant radius
            %   airframe, symmetric fins, and bottom of fins are flush with
            %   bottom of rocket.
            arguments
                % Required
                mass_dry
                XCG_dry
                AF1_length
                AF1_diam
                FN1_root
                FN1_tip
                FN1_height
                FN1_sweep
                THICKNESS
                LMAXU
                NC_length

                % Optional
                NC_diam = AF1_diam
                ZUPPER = (THICKNESS / 2) ./ FN1_root

                % Name-value
                options.LREF = max(NC_diam, AF1_diam)
                options.SREF = pi * (max(NC_diam, AF1_diam) / 2) ^ 2
                options.ROUGH = []
                options.DEXIT = []
                options.FN1_XLE = (AF1_length + NC_length) - FN1_root
                options.ZLOWER = ZUPPER
                options.LMAXL = LMAXU
                options.NC_type = []
                options.NC_power = []
                options.FN1_number = []
                options.FN1_type = []
                options.FN1_NACA = []
                options.FN1_edgeRad = []
                options.m_table = []
                options.COM_table = []
                options.MOI_table = []
                options.MOI_dry = []
                options.Name = []
                options.motor = Motor.empty
                options.main_parachute = Parachute.empty
                options.drogue_parachute = Parachute.empty
            end

            if nargin > 0
                obj.mass_dry = mass_dry;
                obj.XCG_dry = XCG_dry;
                obj.LREF = options.LREF;
                obj.SREF = options.SREF; 
                if ~isempty(options.ROUGH) obj.ROUGH = options.ROUGH; end
                if ~isempty(options.DEXIT) obj.DEXIT = options.DEXIT; end
                if ~isempty(options.NC_type) obj.NC_type = options.NC_type; end
                if obj.NC_type == "POWER" && ~isempty(options.NC_power)
                    obj.NC_power = options.NC_power;
                end
                obj.NC_length = NC_length;
                obj.NC_diam = NC_diam;
                obj.AF1_length = AF1_length;
                obj.AF1_diam = AF1_diam;
                if ~isempty(options.FN1_number) obj.FN1_number = options.FN1_number; end
                if ~isempty(options.FN1_type) obj.FN1_type = options.FN1_type; end
                if isempty(options.FN1_NACA) options.FN1_NACA = obj.FN1_NACA; end
                obj.FN1_NACA(1:numel(obj.FN1_type)) = options.FN1_NACA;
                for f = 1:numel(obj.FN1_type)
                    if obj.FN1_type(f) == "NACA" && ~isempty(options.FN1_NACA(f))
                        obj.FN1_NACA(f) = options.FN1_NACA(f);
                    end
                end
                obj.FN1_root = FN1_root;
                obj.FN1_tip = FN1_tip;
                obj.FN1_height = FN1_height;
                obj.FN1_sweep = FN1_sweep;
                if isempty(options.FN1_edgeRad) options.FN1_edgeRad = obj.FN1_edgeRad; end
                obj.FN1_edgeRad(1:numel(obj.FN1_type)) = options.FN1_edgeRad;
                obj.FN1_XLE = options.FN1_XLE;
                obj.THICKNESS = THICKNESS;
                obj.ZUPPER = ZUPPER;
                obj.ZLOWER = options.ZLOWER;
                obj.LMAXU = LMAXU;
                obj.LMAXL = options.LMAXL;
                obj.LFLATU = 1 - 2 * obj.LMAXU;
                obj.LFLATL = 1 - 2 * obj.LMAXL;
                if ~isempty(options.m_table) obj.m_table = options.m_table; end
                if ~isempty(options.COM_table) obj.COM_table = options.COM_table; end
                if ~isempty(options.MOI_table) obj.MOI_table = options.MOI_table; end
                if ~isempty(options.MOI_dry) obj.MOI_dry = options.MOI_dry; end
                if ~isempty(options.Name) obj.Name = options.Name; end
                if ~isempty(options.motor) obj.motor = options.motor; end
                if ~isempty(options.main_parachute) obj.main_parachute = options.main_parachute; end
                if ~isempty(options.drogue_parachute) obj.drogue_parachute = options.drogue_parachute; end
            end
        end

        function Tx = get.Tx(obj)            
            Tx = obj.motor.thrust_table;
        end

        function obj = modify(obj, modified_properties)
            %MODIFY Modify properties of this class. Must be a table.

            for p = 1:size(modified_properties, 2)
                name = string(modified_properties.Properties.VariableNames(p));
                
                % Perform checks and make corrections if necessary
                % Fix the distance from the nozzle to the leading edge of
                % the fin when changing NC_length so fins don't go off rocket
                if name == "NC_length"
                    NC_length_change = modified_properties.NC_length - obj.NC_length;
                    obj.FN1_XLE = obj.FN1_XLE + NC_length_change;
                    obj.XCG_dry = obj.XCG_dry + NC_length_change;
                    obj.COM_table(:, 2) = obj.COM_table(:, 2) + NC_length_change;
                end
                
                obj.(name) = modified_properties.(name);
            end
        end

        function obj = remove_finset(obj, fin_set)
            %REMOVE_FINSET Remove one of the finsets from the rocket
            
            i = 1:numel(obj.FN1_number);

            fs = ~ismember(i, fin_set);

            obj.FN1_number = obj.FN1_number(fs);
            obj.FN1_type = obj.FN1_type(fs);
            obj.FN1_NACA = obj.FN1_NACA(fs);
            obj.FN1_root = obj.FN1_root(fs);
            obj.FN1_tip = obj.FN1_tip(fs);
            obj.FN1_height = obj.FN1_height(fs);
            obj.FN1_sweep = obj.FN1_sweep(fs);
            obj.FN1_edgeRad = obj.FN1_edgeRad(fs);
            obj.FN1_XLE = obj.FN1_XLE(fs);
            obj.THICKNESS = obj.THICKNESS(fs);
            obj.ZUPPER = obj.ZUPPER(fs);
            obj.ZLOWER = obj.ZLOWER(fs);
            obj.LMAXU = obj.LMAXU(fs);
            obj.LMAXL = obj.LMAXL(fs);
            obj.LFLATU = obj.LFLATU(fs);
            obj.LFLATL = obj.LFLATL(fs);
    end

        function obj = approximate_COM_dry(obj)
            %APPROXIMATE_DRY_MOI Calculate a rough approximation of the dry
            % moment of inertia if there is no data from CAD or experiment.
            
            % ADD IN FINS - seperate class?
            obj.XCG_dry = (obj.NC_length + obj.AF1_length) / 2;
        end

        function obj = approximate_MOI_dry(obj)
            %APPROXIMATE_DRY_MOI Calculate a rough approximation of the dry
            % moment of inertia if there is no data from CAD or experiment.

            length_rocket = (obj.NC_length + obj.AF1_length) * 0.0254; % [m]
            
            % ADD IN FINS - seperate class?
            % Assuming Symmetric Inertia Matrix 
            Ixx = 1 / 2 * obj.mass_dry * (obj.LREF / 2 * 0.0254) ^ 2;
            Ixy = 0;
            Ixz = 0;
            Iyy =  1 / 12 * obj.mass_dry * length_rocket ^ 2; 
            Iyz = 0;
            Izz =  1 / 12 * obj.mass_dry * length_rocket ^ 2; 
            
            obj.MOI_dry = [Ixx, 0, 0;
                           0, Iyy, 0;
                           0, 0, Izz];
        end

        function obj = calculate_COM_and_MOI(obj) % - need to move motor calcs to motor class
            %CALCULATE_COM_and_MOI Calculate center of mass and moment of inertia timeseries from thrust and mass tables properties along with the motor object.
            
            if zeros(3) == obj.XCG_dry
                obj.approximate_COM_dry()
            end

            if isempty(obj.MOI_dry)
                obj.approximate_MOI_dry()
            end

            dt = obj.Tx(2:end, 1) - obj.Tx(1:(end - 1), 1);

            thrust = obj.Tx(2:end, 2);

            impulse = thrust .* dt;
            
            %% Calculate Propellant Mass
            mass_flow = impulse ./ obj.motor.total_impulse .* obj.motor.total_prop_mass ./ dt;

            mass_burn = mass_flow .* dt;

            mass_burn_p975 = 0.975 * mass_burn; % Why?

            prop_mass = obj.motor.total_prop_mass - cumsum(mass_burn_p975);

            mass = max(prop_mass, 0) + obj.motor.dry_mass;

            if isempty(obj.m_table)
                rocket_mass = mass + obj.mass_dry;
                obj.m_table = [obj.Tx(:, 1), [rocket_mass; rocket_mass(end)]];
            end

            %% Center of Mass
            % Structure mass * radius
            bulkhead_mr = obj.motor.bulkhead_m * obj.motor.bulkhead_r;
            casing_mr = obj.motor.casing_m * obj.motor.casing_r;
            liner_mr = obj.motor.liner_m * obj.motor.liner_r;
            nozzle_mr = obj.motor.nozzle_m * obj.motor.nozzle_r;

            dry_mr = bulkhead_mr + casing_mr + liner_mr + nozzle_mr;

            prop_length = (obj.motor.nozzle_r - obj.motor.bulkhead_r) ...
                          .* (1 - cumsum(impulse) ./ sum(impulse) .* 0);

            r_prop = prop_length / 2 + obj.motor.bulkhead_r;

            prop_mr = max(prop_mass .* r_prop, 0);

            COM_prop = (dry_mr + prop_mr) ./ mass;

            COM = (obj.XCG_dry * 0.0254 .* obj.mass_dry + ...
                ((obj.NC_length + obj.AF1_length) * 0.0254 - (obj.motor.motor_length - COM_prop)) .* mass) ...
                  ./ (obj.mass_dry + mass);

            obj.COM_table = [obj.Tx(:, 1), [COM; COM(end)]];

            %% Moment of Inertia
            Ixx_prop = obj.m_table(:, 2) .* obj.motor.prop_radius .^ 2 / 2;
            Iyy_prop = obj.m_table(:, 2) .* (obj.motor.prop_radius .^ 2 / 4 ...
                       + [prop_length; prop_length(end)] .^ 2 / 12 + ([COM_prop; COM_prop(end)] ...
                       + obj.motor.bulkhead_r + [prop_length; prop_length(end)] / 2) .^ 2);
            Izz_prop = Iyy_prop;

            zero = zeros(size(Ixx_prop));

            I_prop = reshape([Ixx_prop, zero, zero; ...
                             zero, Iyy_prop, zero; ...
                             zero, zero, Izz_prop], [], 3, 3);
            
            axis_shift = ([mass; mass(end)] .* (((obj.NC_length + obj.AF1_length) * 0.0254 ...
                            - [COM_prop; COM_prop(end)]) - obj.XCG_dry * 0.0254) .^ 2);

            obj.MOI_table = I_prop ...
                            + reshape([zero, zero, zero; ...
                                       zero, axis_shift, zero; ...
                                       zero, zero, axis_shift], [], 3, 3) ...
                            + reshape(obj.MOI_dry, 1, 3, 3);
        end

        function save(obj)
            %SAVE Save rocket to 'Rockets' folder in a convenient way.
            
            rocket = obj;
            
            path = string("Rockets\" + obj.Name + ".mat");

            save(path, "rocket")
        end

        function [] = display_geo(obj)
            %DISP Display wireframe of rocket geometry
            %   Show a representation of the rocket to verify that it looks
            %   practical and valid.
            
            fin_set_number = numel(obj.FN1_height);
            
            figure
            subplot(1 + 2 * fin_set_number,1,1)
            NC = polyshape([0,     1,      1] *  obj.NC_length, ...
                           [0, 1 / 2, -1 / 2] *  obj.NC_diam); % Need to actually use correct nosecone shape
            AF1 = polyshape([0,                          1,    1,    0] *  obj.AF1_length +  obj.NC_length, ...
                            [ obj.NC_diam / 2,  obj.AF1_diam / 2, - obj.AF1_diam / 2, - obj.NC_diam / 2]);
            for f = 1:fin_set_number
                FN1_top(f) = polyshape([0,  obj.FN1_height(f) * tand(obj.FN1_sweep(f)),  obj.FN1_height(f) * tand(obj.FN1_sweep(f)) +  obj.FN1_tip(f),  obj.FN1_root(f)] +  obj.FN1_XLE(f), ...
                                    [0, 1, 1, 0] *  obj.FN1_height(f) +  obj.AF1_diam / 2);
                if obj.FN1_number(f) > 1
                    FN1_bottom(f) = polyshape([0,  obj.FN1_height(f) * tand(obj.FN1_sweep(f)),  obj.FN1_height(f) * tand(obj.FN1_sweep(f)) +  obj.FN1_tip(f),  obj.FN1_root(f)] +  obj.FN1_XLE(f), ...
                                       -([0, 1, 1, 0] *  obj.FN1_height(f) +  obj.AF1_diam / 2));
                end
            end
            plot(NC); hold on
            plot(AF1); hold on
            p = 2 + fin_set_number;
            for f = 1:fin_set_number
                plot(FN1_top(f)); hold on
                if obj.FN1_number(f) > 1
                    plot(FN1_bottom(f)); hold on
                    p = p + 1;
                end
            end
            scatter(obj.COM_table(:, 2), zeros([size(obj.COM_table, 1), 1]), 15, "blue", "filled"); hold on
            scatter(obj.XCG_dry, 0,30,"red", "filled"); hold off
            axis equal
            xlabel("Axial [m]")
            ylabel("Transverse [m]")
            title("Rocket Diagram")
            legend([strings(1, p), "COM table", "XCG Dry"])
            
            for f = 1:fin_set_number
                subplot(1 + 2 * fin_set_number,1,2 * f)
                if obj.FN1_type(f) == "HEX"
                    ZUPPERs = [obj.ZUPPER(f), obj.ZUPPER(f) * obj.FN1_root(f) / obj.FN1_tip(f)];
                    ZLOWERs = [obj.ZLOWER(f), obj.ZLOWER(f) * obj.FN1_root(f) / obj.FN1_tip(f)];
                    LMAXUs = [obj.LMAXU(f), obj.LMAXU(f) * obj.FN1_root(f) / obj.FN1_tip(f)];
                    LMAXLs = [obj.LMAXL(f), obj.LMAXL(f) * obj.FN1_root(f) / obj.FN1_tip(f)];
                    LFLATUs = [obj.LFLATU(f),  1 - 2 * LMAXUs(2)];
                    LFLATLs = [obj.LFLATL(f),  1 - 2 * LMAXLs(2)];
    
                    FN1_root_shape = polyshape([0,  LMAXUs(1),  LMAXUs(1) +  LFLATUs(1), 1,  LMAXLs(1) +  LFLATLs(1),  LMAXLs(1)] *  obj.FN1_root(f), ...
                                        [0,  ZUPPERs(1),  ZUPPERs(1), 0, - ZLOWERs(1), - ZLOWERs(1)] *  obj.FN1_root(f));
                elseif obj.FN1_type(f) == "ARC"
                    points = 15;
                    x = linspace(0, 1, points + 2);
                    x = x(2:(end - 1));
                    y_upper = obj.ZUPPER(f) * sin(pi * x);
                    y_lower = -obj.ZLOWER(f) * sin(pi * x);
                    FN1_root_shape = polyshape([0, x, 1, flip(x)] *  obj.FN1_root(f), ...
                                        [0, y_upper, 0, y_lower] *  obj.FN1_root(f));
                end
                plot(FN1_root_shape);
                axis equal
                xlabel("Axial [m]")
                ylabel("Transverse [m]")
                title(sprintf("Fin Set %g Cross Section at Root", f))
                
                subplot(1 + 2 * fin_set_number,1,2 * f + 1)
                if obj.FN1_type(f) == "HEX"
                    FN1_tip_shape = polyshape([0,  LMAXUs(2),  LMAXUs(2) +  LFLATUs(2), 1,  LMAXLs(2) +  LFLATLs(2),  LMAXLs(2)] *  obj.FN1_tip(f), ...
                                        [0,  ZUPPERs(2),  ZUPPERs(2), 0, - ZLOWERs(2), - ZLOWERs(2)] *  obj.FN1_tip(f));
                elseif obj.FN1_type(f) == "ARC"
                    points = 15;
                    x = linspace(0, 1, points + 2);
                    x = x(2:(end - 1));
                    y_upper = obj.ZUPPER(f) * sin(pi * x);
                    y_lower = -obj.ZLOWER(f) * sin(pi * x);
                    FN1_tip_shape = polyshape([0, x, 1, flip(x)] *  obj.FN1_tip(f), ...
                                        [0, y_upper, 0, y_lower] *  obj.FN1_tip(f));
                end
                plot(FN1_tip_shape);
                axis equal
                xlabel("Axial [m]")
                ylabel("Transverse [m]")
                title(sprintf("Fin Set %g Cross Section at Tip", f))
            
            end

            sgtitle("Geometry Plots for " + obj.Name)
            hold off
        end

        function aerotable = create_aerotable(rocket, thread_number, target_batch_size, options)
            %CREATE_AEROTABLE Create aerotable for the rocket and save
            %   summary
            arguments
                rocket
                thread_number

                target_batch_size = 100

                options.ALPHA = [0,0.25,0.5,.75,1,1.25,1.5,1.75,2,2.25,2.5,3,3.5,4,4.5,5,6,7,10,20] % Array of angles of attack, (deg)
                options.MACH = [.01,0.1,0.3,0.7,0.9,.95,.98,1.01,1.02,1.04,1.07,1.1,1.2,1.3,1.5,1.75,2.0,2.5,3.0,3.5] % Array of Mach numbers
                options.ALT_iterations = 40 % Number of altitude breakpoints
                options.ALT_max = 50000 % Max rocket altitude, (ft)
                options.PHI_iterations = 40 % Number of phi breakpoints
                options.aerotable_name = ""
                options.Name_tag = ""
                options.save_aerotable = true
                options.aerotable_parameter_dialog = false;
            end
            
            if options.aerotable_parameter_dialog
                prompt = {'Enter Max ALT [ft]:', 'Enter ALT Iterations:', 'Enter PHI Iterations', 'Enter ALPHA [deg]', 'Enter MACH'};
                dlgtitle = 'Aerotable Parameters';
                fieldsize = [1 45; 1 45; 1 45; 1 70; 1 70];
                definput = {char(string(options.ALT_max)), char(string(options.ALT_iterations)), char(string(options.PHI_iterations)), mat2str(options.ALPHA), mat2str(options.MACH)};
                answer = inputdlg(prompt,dlgtitle,fieldsize,definput);

                if ~isempty(answer)
                    options.ALT_max = str2double(answer{1});
                    options.ALT_iterations = str2double(answer{2});
                    options.PHI_iterations = str2double(answer{3});
                    options.ALPHA = eval(answer{4});
                    options.MACH = eval(answer{5});
                else
                    error("Exited Aerotable Creation");
                end
            end

            aerotable_parameters = struct();
    
            % Specify the structure of the aerotables
            % DATCOM expects MACH and ALPHA to have maximum sizes of 20 
            aerotable_parameters.ALPHA = options.ALPHA;
            aerotable_parameters.MACH = options.MACH;
            aerotable_parameters.ALT_iterations = options.ALT_iterations;
            aerotable_parameters.ALT_max = options.ALT_max;
            aerotable_parameters.PHI_iterations = options.PHI_iterations;

            aerotable = generate_aerotables(rocket, aerotable_parameters, target_batch_size = target_batch_size, thread_number = thread_number);
        
            
            if options.aerotable_name == ""
                aerotable_name = sprintf("%s_Aerotable_%s-%gALT_%gPHI", rocket.Name, string(options.ALT_max), options.ALT_iterations, options.PHI_iterations);
            else
                aerotable_name = options.aerotable_name;
            end

            aerotable_name = aerotable_name + options.Name_tag;

            while options.save_aerotable && isfile(string("Aero Databases\" + aerotable_name + ".mat"))
                answer = questdlg(sprintf('Aerotable with name %s already exists. What should be done?', aerotable_name), ...
                    'Aerotable Save Conflict', ...
                    "Replace","Don't save", "Append Tag", "Append Tag");
                
                switch answer
                    case "Replace"
                        break;
                    case "Don't save"
                        options.save_aerotable = false;
                    case "Append Tag"
                        prompt = {sprintf('Enter Tag for %s:', aerotable_name)};
                        dlgtitle = 'Aerotable Tag';
                        fieldsize = [1 45];
                        definput = {'_2'};
                        answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
                        
                        if ~isempty(answer)
                            aerotable_name = aerotable_name + string(answer);
                        else
                            options.save_aerotable = false;
                        end
                end
            end

            if options.save_aerotable
                save(string("Aero Databases\" + aerotable_name), "aerotable")
            end
        end
        
        function r = plus(upperstage, lowerstage)
            rocket_COM_table = (upperstage.COM_table(1, 2) .* upperstage.m_table(1, 2) + (lowerstage.COM_table + upperstage.NC_length + upperstage.AF1_length) .* lowerstage.m_table) ./ (upperstage.m_table(1, 2) .* lowerstage.m_table);
            rocket_COM_table(:, 1) = lowerstage.COM_table(:, 1);
            rocket_COM_dry = (upperstage.XCG_dry .* upperstage.mass_dry + (lowerstage.XCG_dry + upperstage.NC_length + upperstage.AF1_length) .* lowerstage.mass_dry) ./ (upperstage.mass_dry .* lowerstage.mass_dry);

            upperstage_COM_dist = interp1(rocket_COM_table(:, 1), rocket_COM_table(:, 2) - upperstage.COM_table(1, 2), linspace(rocket_COM_table(1, 1), rocket_COM_table(end, 1), numel(lowerstage.MOI_table(:, 1, 1))))';
            lowerstage_COM_dist = interp1(rocket_COM_table(:, 1), (lowerstage.COM_table(:, 2) + upperstage.NC_length + upperstage.AF1_length) - rocket_COM_table(:, 2), linspace(rocket_COM_table(1, 1), rocket_COM_table(end, 1), numel(lowerstage.MOI_table(:, 1, 1))))';

            resampled_lowerstage_mass = interp1(lowerstage.m_table(:, 1), lowerstage.m_table(:, 2), linspace(lowerstage.m_table(1, 1), lowerstage.m_table(end, 1), numel(lowerstage.MOI_table(:, 1, 1))))';

            rocket_MOI_table = upperstage.MOI_table(1, :, :) + lowerstage.MOI_table;
            rocket_MOI_table(:, 2, 2) = rocket_MOI_table(:, 2, 2) + upperstage.m_table(1, 2) .* upperstage_COM_dist .^ 2 + resampled_lowerstage_mass .* lowerstage_COM_dist .^ 2;
            rocket_MOI_table(:, 3, 3) = rocket_MOI_table(:, 3, 3) + upperstage.m_table(1, 2) .* upperstage_COM_dist .^ 2 + resampled_lowerstage_mass .* lowerstage_COM_dist .^ 2;
    
            upperstage_COM_dry_dist = rocket_COM_dry - upperstage.XCG_dry;
            lowerstage_COM_dry_dist = (lowerstage.XCG_dry + upperstage.NC_length + upperstage.AF1_length) - rocket_COM_dry;

            rocket_MOI_dry = upperstage.MOI_dry + lowerstage.MOI_dry;
            rocket_MOI_dry(2, 2) = rocket_MOI_dry(2, 2) + upperstage.mass_dry .* upperstage_COM_dry_dist .^ 2 + lowerstage.mass_dry .* lowerstage_COM_dry_dist .^ 2;
            rocket_MOI_dry(3, 3) = rocket_MOI_dry(3, 3) + upperstage.mass_dry .* upperstage_COM_dry_dist .^ 2 + lowerstage.mass_dry .* lowerstage_COM_dry_dist .^ 2;

            r = Rocket(upperstage.mass_dry + lowerstage.mass_dry, ...
                       rocket_COM_dry, ...
                       upperstage.AF1_length + lowerstage.AF1_length, ...
                       lowerstage.AF1_diam, ...
                       [upperstage.FN1_root, lowerstage.FN1_root], ...
                       [upperstage.FN1_tip, lowerstage.FN1_tip], ...
                       [upperstage.FN1_height, lowerstage.FN1_height], ...
                       [upperstage.FN1_sweep, lowerstage.FN1_sweep], ...
                       [upperstage.THICKNESS, lowerstage.THICKNESS], ...
                       [upperstage.LMAXU, lowerstage.LMAXU], ...
                       upperstage.NC_length, ...
                       FN1_XLE = [upperstage.FN1_XLE, lowerstage.FN1_XLE + upperstage.NC_length - lowerstage.NC_length + upperstage.AF1_length], ...
                       FN1_number = [upperstage.FN1_number, lowerstage.FN1_number], ...
                       FN1_type = [upperstage.FN1_type, lowerstage.FN1_type], ...
                       NC_type = upperstage.NC_type, ...
                       Name = sprintf("%s%s", upperstage.Name, lowerstage.Name), ...
                       m_table = [0, upperstage.m_table(1, 2)] + lowerstage.m_table, ...
                       COM_table = rocket_COM_table, ...
                       MOI_table = rocket_MOI_table, ...
                       MOI_dry = rocket_MOI_dry, ...
                       motor = lowerstage.motor);
        end
    end
end

