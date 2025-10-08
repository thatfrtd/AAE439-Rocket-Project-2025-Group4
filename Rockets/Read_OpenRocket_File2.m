function [parsedOpenRocket, launch_site] = Read_OpenRocket_File2(openrocket_file, options)
%READ_OPENROCKET_FILE2 Reads an OpenRocket file as an XML file and creates a motor and rocket classes to be used in the 6DoF
%   The file must contain a simulation with the desired motor being used in
%   the primary simulation (sim corresponding to primary_sim_number). 
%   openrocket_file must contain '.ork' extension.
%
% TODO: 
% - Add more options for importing different fin cross sections
% - Import recovery information
% - Import staging and ignition information
% - Support transitions
% - Import protruberances
% - Support motor tubes extending out of rocket - need to make it only
% affect lengths of individual stage 
% - Use this way of importing OR data for adding noise to OR sims
% - Import multiple different motor configurations
arguments
    openrocket_file

    options.read_tables = true
    options.show_xml = false
    options.plot_rocket_geometry = true
    options.LMAXU = 0.05
    options.rocket_name = ""
    options.primary_sim_number = 1
    options.launch_site_name = []
    options.single_fins_are_protruberances = true
    options.generate_aerotable = true
end

%% Import OpenRocket Data
% Unzip and read as xml
full_openrocket_file = which(openrocket_file);
[filepath, openrocket_file_name, ext] = fileparts(full_openrocket_file);
status = copyfile(full_openrocket_file, openrocket_file_name + ".zip");
unzipped_ork = unzip(fullfile(pwd, openrocket_file_name + ".zip"));
docNode = xmlread(unzipped_ork{1});

% Display xml
if options.show_xml
    xmlwrite(docNode)
end

% Delete zip folder
delete(openrocket_file_name + ".zip")

r = readstruct("rocket.ork",FileType="xml");

stage_number = numel(r.rocket.subcomponents.stage);

% Rocket Name
if options.rocket_name == ""
    rocketName = r.rocket.name;
else
    rocketName = options.rocket_name;
end

for stage_num = 1:stage_number
    s = r.rocket.subcomponents.stage(stage_num);
    
    stage_name = s.name;
    
    use_overrides = false;

    stage_fields = string(fieldnames(s));
    % Rocket Mass
    if ismember("overridemass", stage_fields) && use_overrides
        rocketMass = s.overridemass;
    else
        rocketMass = -1;
    end
    
    % Rocket CG
    if ismember("overridecg", stage_fields) && use_overrides
        rocketCG = s.overridecg;
    else
        rocketCG = -1;
    end
    
    stage_subcomponents_fields = string(fieldnames(s.subcomponents));
    % Nosecone
    if ismember("nosecone", stage_subcomponents_fields)
        nosecone = s.subcomponents.nosecone;
        rocketNCLength = nosecone.length;

        % Rocket NC type
        oprock2datcom_NoseType = struct();
        oprock2datcom_NoseType.karman = "KARMAN";
        oprock2datcom_NoseType.haak = "HAACK";
        oprock2datcom_NoseType.ogive = "OGIVE";
        oprock2datcom_NoseType.conical = "CONICAL";
        oprock2datcom_NoseType.power = "POWER";
        
        rocketNCopType = nosecone.shape;
        
        rocketNCShapeParam = nosecone.shapeparameter;
        
        if rocketNCopType == "haack"
            % Assuming Karman (shape parameter = 0) or L-V Haak (shape parameter = 0.333) noscone
            % Would have to make custom nosecone
            rocketNCShapeParam = round(rocketNCShapeParam * 3);
        
            if rocketNCShapeParam == 0
                rocketNCType = "KARMAN";
            elseif rocketNCShapeParam == 1
                rocketNCType = "HAAK";
            end
        else
            % Assuming Tangent Ogive (shape parameter == 1) if ogive. Would have to
            % make custom nosecone
            rocketNCType = oprock2datcom_NoseType.(rocketNCopType);
        end
    else
        rocketNCLength = 0;
        rocketNCType = [];
    end
    
    bodytube = s.subcomponents.bodytube;
    bodytube_diam = 2 * remove_auto([bodytube.radius]);
    bodytube_length = remove_auto([bodytube.length]);
    
    % Assume no transitions in body - ADD IN ASAP
    rocketAFLength = sum(bodytube_length);
    rocketAFDiameter = bodytube_diam(1);
    
    % Fins
    fin_number = [];
    fin_root = [];
    fin_tip = [];
    fin_height = [];
    fin_sweep = [];
    fin_thickness = [];
    fin_xle = [];
    
    rocketAFLengthf = 0;
    for l = 1:numel(bodytube)
        rocketAFLength_l = bodytube(l).length;
    
        if ismember("trapezoidfinset", string(fieldnames(bodytube(l).subcomponents)))
            trapfinset_number = numel(bodytube(l).subcomponents.trapezoidfinset);
            trapnodes = bodytube(l).subcomponents.trapezoidfinset;

            if options.single_fins_are_protruberances
                fincounts = [trapnodes.fincount];

                trapfinset_number = trapfinset_number - sum(fincounts == 1);
                trapnodes = trapnodes(fincounts ~= 1);

                % Add protruberance info somehow - TODO
            end
        else
            trapfinset_number = 0;
        end
        
        trap_fn_number = zeros([1, trapfinset_number]);
        trap_fn_root = zeros([1, trapfinset_number]);
        trap_fn_tip = zeros([1, trapfinset_number]);
        trap_fn_height = zeros([1, trapfinset_number]);
        trap_fn_sweep = zeros([1, trapfinset_number]);
        trap_fn_xle = zeros([1, trapfinset_number]);
        trap_fn_thickness = zeros([1, trapfinset_number]);
        
        for tf = 1:trapfinset_number
            trapnode = trapnodes(tf);

            % Rocket FN1 number
            trap_fn_number(tf) = trapnode.fincount;
            
            trap_fn_root(tf) = trapnode.rootchord;
            trap_fn_tip(tf) = trapnode.tipchord;
            trap_fn_height(tf) = trapnode.height;
            trap_fn_sweep_length = trapnode.sweeplength;
            trap_fn_sweep(tf) = atan2d(trap_fn_sweep_length, trap_fn_height(tf));
            
            % Fin position type
            rocketFNpositiontype = trapnode.position.typeAttribute;
            
            % Rocket FN1 XLE
            trap_fn_xle(tf) = trapnode.position.Text;
            
            if rocketFNpositiontype == "bottom"
                trap_fn_xle(tf) = trap_fn_xle(tf) - trap_fn_root(tf) + rocketAFLengthf + rocketAFLength_l + rocketNCLength;
            elseif rocketFNpositiontype == "top"
                trap_fn_xle(tf) = trap_fn_xle(tf) - trap_fn_root(tf) + rocketAFLengthf + rocketNCLength;
            elseif rocketFNpositiontype == "middle"
                trap_fn_xle(tf) = trap_fn_xle(tf) - trap_fn_root(tf) + rocketAFLengthf + rocketAFLength_l / 2 + rocketNCLength;
            end
            
            % Rocket fin thickness
            trap_fn_thickness(tf) = trapnode.thickness;
        end
        
        if ismember("freeformfinset", string(fieldnames(bodytube(l).subcomponents)))
            freeformfinset_number = numel(bodytube(l).subcomponents.freeformfinset);
            freeformnodes = bodytube(l).subcomponents.freeformfinset;

            if options.single_fins_are_protruberances
                fincounts = [freeformnodes.fincount];

                freeformfinset_number = freeformfinset_number - sum(fincounts == 1);
                freeformnodes = freeformnodes(fincounts ~= 1);

                % Add protruberance info somehow - TODO
            end
        else
            freeformfinset_number = 0;
        end
    
        free_fn_number = zeros([1, freeformfinset_number]);
        free_fn_root = zeros([1, freeformfinset_number]);
        free_fn_tip = zeros([1, freeformfinset_number]);
        free_fn_height = zeros([1, freeformfinset_number]);
        free_fn_sweep = zeros([1, freeformfinset_number]);
        free_fn_xle = zeros([1, freeformfinset_number]);
        free_fn_thickness = zeros([1, freeformfinset_number]);
        
        for ff = 1:freeformfinset_number
            freeformnode = freeformnodes(ff);

            % Rocket FN1 number
            free_fn_number(ff) = freeformnode(ff).fincount;
            
            % Rocket FN1 root, tip, height, sweep
            % Approximate the freeform fin in a way DATCOM likes
            rocketFNPointsx = [freeformnode.finpoints.point.xAttribute];
            rocketFNPointsy = [freeformnode.finpoints.point.yAttribute];
            
            % Assume there are 4 points making up the fin
            finpoints = zeros([4, 2]);
            
            for p = 1:4
                finpoints(p, 1) = rocketFNPointsx(p);
                finpoints(p, 2) = rocketFNPointsy(p);
            end
            
            free_fn_root(ff) = finpoints(4, 1);
            free_fn_tip(ff) = finpoints(3, 1) - finpoints(2,1);
            free_fn_height(ff) = sum(finpoints(2:3, 2)) / 2;
            free_fn_sweep(ff) = atan2d(finpoints(2,1), finpoints(2,2));
            
            % Fin position type
            rocketFNpositiontype = freeformnode.position.typeAttribute;
            
            % Rocket FN1 XLE
            free_fn_xle(ff) = freeformnode.position.Text;
            
            if rocketFNpositiontype == "bottom"
                free_fn_xle(ff) = free_fn_xle(ff) - free_fn_root(ff) + rocketAFLengthf + rocketAFLength_l + rocketNCLength;
            elseif rocketFNpositiontype == "top"
                free_fn_xle(ff) = free_fn_xle(ff) + rocketAFLengthf + rocketNCLength;
            end
            
            % Rocket fin thickness
            free_fn_thickness(ff) = freeformnode.thickness;
        end
    
        fin_number = [fin_number, free_fn_number, trap_fn_number];
        fin_root = [fin_root, free_fn_root, trap_fn_root];
        fin_tip = [fin_tip, free_fn_tip, trap_fn_tip];
        fin_height = [fin_height, free_fn_height, trap_fn_height];
        fin_sweep = [fin_sweep, free_fn_sweep, trap_fn_sweep];
        fin_thickness = [fin_thickness, free_fn_thickness, trap_fn_thickness];
        fin_xle = [fin_xle, free_fn_xle, trap_fn_xle];
    
        rocketAFLengthf = rocketAFLengthf + rocketAFLength_l;
    end
    
    [fin_xle, finset_indices] = sort(fin_xle, 2, "ascend");
    
    fin_number = fin_number(finset_indices);
    fin_root = fin_root(finset_indices);
    fin_tip = fin_tip(finset_indices);
    fin_height = fin_height(finset_indices);
    fin_sweep = fin_sweep(finset_indices);
    fin_thickness = fin_thickness(finset_indices);
    fin_type(1:numel(fin_number)) = "HEX";
    LMAXU(1:numel(fin_number)) = options.LMAXU;
    
    transition_array = zeros([numel(uniquetol(bodytube_diam)) + 1, 1]);
    transition_end_diam = [0, 0];

    % Transition
    if ismember("transition", stage_subcomponents_fields)
        transition = s.subcomponents.transition;
        transitionLength = remove_auto([transition.length]);
        transitionFore = 2 * remove_auto([transition.foreradius]);
        transitionAft = 2 * remove_auto([transition.aftradius]);

        transition_i = 1;
        for b = 1:numel(transition_array)
            if transition_i <= numel(transitionLength) && ismembertol(transitionAft(transition_i), bodytube_diam(b))
                transition_array(b) = transitionLength(transition_i);

                if b == 1
                    transition_end_diam(1) = transitionFore(1);
                end

                transition_i = transition_i + 1;
            elseif transition_i <= numel(transitionLength) && ismembertol(transitionFore(transition_i), bodytube_diam(b))
                transition_array(b) = transitionLength(transition_i);
               
                if b == numel(transition_array)
                    transition_end_diam(2) = transitionAft(b);
                end

                transition_i = transition_i + 1;
            end
        end
    end

    % Read in simulation data to get thrust, center of mass, moment of
    % intertia, and mass throughout the burn
    flight = r.simulations.simulation(options.primary_sim_number).flightdata.databranch(stage_num).datapoint;
    col_names = split(r.simulations.simulation(options.primary_sim_number).flightdata.databranch(stage_num).typesAttribute, ",");
    data = zeros([numel(flight), 54]);
    for t = 1:numel(flight) 
        data(t, :) = str2num(flight(t));
    end
    flightdata = array2table(data, "VariableNames",col_names);
    
    flightevents = r.simulations.simulation(options.primary_sim_number).flightdata.databranch(stage_num).event;
    
    ignition_i = 1;
    burnout_i = 1;
    separation_i = 1;
    for e = 1:numel(flightevents)
        if flightevents(e).typeAttribute == "ignition"
            if ignition_i == (stage_number - stage_num + 1)
                t_ignition = flightevents(e).timeAttribute;
            end
            ignition_i = ignition_i + 1;
        end
        
        if flightevents(e).typeAttribute == "burnout"
            if burnout_i == (stage_number - stage_num + 1)
                t_burnout = flightevents(e).timeAttribute;
            end
            burnout_i = burnout_i + 1;
        end

        if flightevents(e).typeAttribute == "stageseparation"
            if separation_i == (stage_number - stage_num + 1)
                t_separation = flightevents(e).timeAttribute;
            end
            separation_i = separation_i + 1;
        end
    end

    burn_period = flightdata.Time >= t_ignition & flightdata.Time <= t_burnout;

    flightdata_burn = flightdata(burn_period, :);
    flightdata_burn.Time = flightdata_burn.Time - t_ignition;

    launchconditions = r.simulations.simulation(options.primary_sim_number).conditions;
    
    if isempty(options.launch_site_name)
        options.launch_site_name = sprintf("Lat:%.2f Long:%.2f", launchconditions.launchlatitude, launchconditions.launchlongitude);
    end

    launch_site = LaunchSite(options.launch_site_name, ...
                             launchconditions.launchlatitude, ...
                             launchconditions.launchlongitude, ...
                             launchconditions.launchrodangle, ...
                             launchconditions.launchroddirection, ...
                             90 - launchconditions.launchaltitude, ...
                             launchconditions.launchrodlength);
    
    % Can't find wind direction in launchconditions...
    %wind = ConstantWind("ORConstantWind", , launchconditions.windaverage)
    
    rocket_mass_table = [flightdata_burn.Time, flightdata_burn.Mass];
    rocket_COM_table = [flightdata_burn.Time, flightdata_burn.("CG location")];
    rocket_XCG_dry = rocket_COM_table(end, 2);
    
    if rocketCG == -1
        rocketCG = rocket_XCG_dry;
    end

    zeroes = zeros([1, numel(flightdata_burn.Time)]);
    rocket_MOI_table = reshape([flightdata_burn.("Rotational moment of inertia")', zeroes, zeroes, ...
                        zeroes, flightdata_burn.("Longitudinal moment of inertia")', zeroes, ...
                        zeroes, zeroes, flightdata_burn.("Longitudinal moment of inertia")'], [numel(flightdata_burn.Time), 3, 3]);
    rocket_MOI_dry = squeeze(rocket_MOI_table(end, :, :));
    
    thrust_table = [flightdata_burn.Time, flightdata_burn.Thrust];
    motor_mass_table = [flightdata_burn.Time, flightdata_burn.("Motor mass")];
    
    motor_dry_mass = motor_mass_table(end, 2);
    propellant_mass = motor_mass_table(1, 2) - motor_mass_table(end, 2);
    
    innertubes = [r.rocket.subcomponents.stage(stage_num).subcomponents.bodytube(end).subcomponents.innertube];
    motor_str = innertubes(1).motormount.motor;
    motor_name = sprintf("%s %s", motor_str(1).manufacturer, motor_str(1).designation);
    motor_length = motor_str(1).length;
    motor_diameter = motor_str(1).diameter;
    
    % Assume only two stage rocket
    % Store the booster alone's COM, mass, MOI
    if stage_num == 2
        separation_i = find(flightdata.Time > t_separation, 1);
        flight_data_booster = flightdata(separation_i:(separation_i + 1), :);

        booster_mass_table = [flight_data_booster.Time, flight_data_booster.Mass];
        booster_COM_table = [flight_data_booster.Time, flight_data_booster.("CG location")] - (openRocket(1).NC_length + openRocket(1).AF1_length);
        booster_XCG_dry = booster_COM_table(end, 2);

        zeroes = zeros([1, numel(flight_data_booster.Time)]);
        booster_MOI_table = reshape([flight_data_booster.("Rotational moment of inertia")', zeroes, zeroes, ...
                            zeroes, flight_data_booster.("Longitudinal moment of inertia")', zeroes, ...
                            zeroes, zeroes, flight_data_booster.("Longitudinal moment of inertia")'], [numel(flight_data_booster.Time), 3, 3]);
        booster_MOI_dry = squeeze(booster_MOI_table(end, :, :));
    end

    %% Construct Motor
    motor = SolidMotor(motor_name, thrust_table, motor_mass_table, [], motor_dry_mass, [], propellant_mass, motor_length, motor_diameter, [], [], [], [], [], []);

    %% Construct Rocket
    openRocket(stage_num) = Rocket(rocketMass, rocketCG, rocketAFLength, rocketAFDiameter, fin_root, fin_tip, fin_height, fin_sweep, fin_thickness, LMAXU, rocketNCLength, FN1_XLE = fin_xle, FN1_number = fin_number, FN1_type = fin_type, NC_type = rocketNCType, Name = stage_name, m_table = rocket_mass_table, COM_table = rocket_COM_table, MOI_table = rocket_MOI_table, MOI_dry = rocket_MOI_dry, motor = motor);
    
    if options.plot_rocket_geometry
        openRocket(stage_num).display_geo();
    end
end

if stage_number > 1
    multistageOpenRocket = MultiStageRocket.stack(rocketName, flip(openRocket));

    % Fix booster and stack... need better way that accounts for it right away
    multistageOpenRocket.stages(1).COM_table = openRocket(2).COM_table;
    multistageOpenRocket.stages(1).MOI_table = openRocket(2).MOI_table;
    multistageOpenRocket.stages(1).MOI_dry = openRocket(2).MOI_dry;
    multistageOpenRocket.stages(1).m_table = openRocket(2).m_table;
    multistageOpenRocket.stages(1).XCG_dry = openRocket(2).XCG_dry;

    multistageOpenRocket.stages(2).COM_table = booster_COM_table;
    multistageOpenRocket.stages(2).MOI_table = booster_MOI_table;
    multistageOpenRocket.stages(2).MOI_dry = booster_MOI_dry;
    multistageOpenRocket.stages(2).m_table = booster_mass_table;
    multistageOpenRocket.stages(2).XCG_dry = booster_XCG_dry;

    if options.plot_rocket_geometry
        multistageOpenRocket.display_geo();
    end

    parsedOpenRocket = multistageOpenRocket;

    if options.generate_aerotable
        aerotables = parsedOpenRocket.create_aerotables(3, 200);
        save(sprintf("Aero Databases/%s_array", parsedOpenRocket.Name), "aerotables");
    end
else
    parsedOpenRocket = openRocket(1);

    if options.generate_aerotable   
        parsedOpenRocket.create_aerotable(3, 200, aerotable_parameter_dialog = true);
    end
end

parsedOpenRocket.save();

% Clean up xml file
delete(unzipped_ork{1});


function [fixed_dimension] = remove_auto(dimension)
    if isstring(dimension)
        fixed_dimension = str2double(erase(dimension, "auto "));
    else
        fixed_dimension = dimension;
    end
end
end