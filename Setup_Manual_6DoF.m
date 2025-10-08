%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HA Sims
% Setup 6DoF Simulink Model for Manual Running Script
% Author: Travis Hastreiter 
% Created On: 3 October, 2023
% Description: Loads in all rocket and aerotable data into the 6DoF
% Simulink model so that it can be manualy run from the Simulink window. If
% the requested aerotable is not found in the Aero Databases folder, it makes 
% a new one. Assigns all variables to the model workspace so that they will 
% persist and not clutter the base workspace.
% Most Recent Change: Travis Hastreiter 14 April, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
model_name = "continuous_6DoF_scoped";

rocket_name = "hitech";
aerotable_name = "hitech_Aerotable_6000-40ALT_40PHI"; % Creates aerotable if it does not exist

launch_site_name = "KSFO";
%wind_name = "ConstantNoWind";
wind_name = "HWM_turb";

%% Launch Parameters
launch_site = load(launch_site_name).launch_site;
launch_site.elevation_angle = 80; % [deg] Angle from horizontal
launch_site.compass_angle = 0; % [deg] Angle from North (+90 degrees is East)

% Wind
wind = load(wind_name).wind;
%wind.start_time = datetime("now");
wind.include_turbulence = false;
%wind.wind_speed = 5.8; % [m/s] Constant wind speed
%wind.wind_direction = 315; % [deg] Direction that wind comes from (+90 degrees is East)

%% Load Rocket
rocket = load(rocket_name).rocket;

%rocket.main_parachute.CD = 0;
rocket.main_parachute.deploy_height = 50;
rocket.drogue_parachute.CD = 0;
% Look at the rocket geometry
%rocket.display_geo()

%% Load Aerotable
if isfile(string("Aero Databases/" + aerotable_name + ".mat"))
    aerotable = load(aerotable_name).aerotable;
else
    warning("Could not find requested aerotable, making new one.")
    
    thread_number = 1;
    aerotable = rocket.create_aerotable(thread_number, 1, aerotable_name = aerotable_name, aerotable_parameter_dialog = true);
end

%% Initialize 6DoF
% Run rail sim
rail_inputs = initialize_rail(rocket, launch_site);
rail_outputs = Run_Simulation("Rail_Sim", rail_inputs);

% Create model inputs
initial_conditions = rail_output_to_6DoF_initial_conditions(rail_outputs, launch_site);
model_inputs = initialize_6DoF(rocket, wind, initial_conditions, aerotable, Apogee_stop=false, Burnout_stop=false);

% Open model so input can be filled and it can be ran
open(string("6DOF Models/" + model_name));

% Send data in model_input struct to simulation
mdlWks = get_param(model_name, 'ModelWorkspace');

input_field_names = fieldnames(model_inputs);
for i = 1:length(input_field_names)
    name = string(input_field_names(i));

    mdlWks.assignin(name,model_inputs.(name));
end
mdlWks.assignin("Rocket", rocket)
mdlWks.assignin("Aerotable", aerotable)