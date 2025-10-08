function [model_inputs] = initialize_6DoF(rocket, wind, initial_conditions, model_inputs, options)
%INTIALIZATION_FUNCTION Initializes 6DoF with provided rocket parameters
%   Store all initialized variables in input struct so it can be feed to
%   the 6DoF.
% initial_conditions:
% time0
% HG0
% x0_NED
% launch_elevation
% uvw0_B
% BG0
% lambda0
% mu0
% w0_BG

arguments
    rocket {mustBeA(rocket, "Rocket")}
    wind {mustBeA(wind, "Wind")}
    initial_conditions
    model_inputs = struct() % Other model inputs can be passed in here (like an aerotable)
    options.Apogee_stop = false
    options.Burnout_stop = false
    options.turbulence_noise_seeds (1, 4) = [23341 23342 23343 23344];
end

%% Physical Parameters
g = 9.81; % [m/s2] Gravitational Acceleration at Sea Level
w_E = 2 * pi / (86164); % [rad/s] Earth's Rotation
R_earth = 6371 * 10 ^ 3; % [m] Earth's Radius

%% Simulink Initilizations
GB0 = initial_conditions.BG0.'; % Body to Geodetic

v0_NED = GB0*initial_conditions.uvw0_B; % Initial velocity in North East Down frame

latlong0 = [initial_conditions.lambda0, initial_conditions.mu0];

% Initial angular velociy in body frame
w0_BN = initial_conditions.w0_BG + initial_conditions.BG0 ...
    * [w_E*cos(initial_conditions.lambda0);
       0;
       -w_E*sin(initial_conditions.lambda0)];

%% Put variables into model inputs struct
% Store model inputs
model_inputs.uvw0_B = initial_conditions.uvw0_B;
model_inputs.v0_NED = v0_NED;
model_inputs.x0_NED = initial_conditions.x0_NED;
model_inputs.Apogee_stop = options.Apogee_stop;
model_inputs.Burnout_stop = options.Burnout_stop;
model_inputs.HG0 = initial_conditions.HG0;
model_inputs.latlong0 = latlong0;
model_inputs.w0_BN = w0_BN;
model_inputs.BG0 = initial_conditions.BG0;
model_inputs.g = g;
model_inputs.R_earth = R_earth;
model_inputs.w_E = w_E;
model_inputs.time0 = initial_conditions.time0;
model_inputs.launch_elevation = initial_conditions.launch_elevation;

if wind.include_turbulence
    % Set the 'wingspan' used to calculate the angular velocities from
    % turbulence
    % Adjust for when there are an odd amount of fins?
    model_inputs.finspan = rocket.LREF + 2 * max(rocket.FN1_sweep);

    % Wind velocity above 6m ground level is used low level turbulence
    model_inputs.wind_speed6m = 5; % REPLACE WITH ACTUAL VALUES - make evaluate abstract method for Wind?
    model_inputs.wind_direction6m = 0; % REPLACE WITH ACTUAL VALUES

    model_inputs.turbulence_noise_seeds = options.turbulence_noise_seeds;
end

% Copy necessary wind data
model_inputs = copy_class_fields_to_struct(wind, model_inputs);
% Copy parachute data
model_inputs = copy_class_fields_to_struct(rocket.main_parachute, model_inputs, "main_");
model_inputs = copy_class_fields_to_struct(rocket.drogue_parachute, model_inputs, "drogue_");

% Give the model the requisite rocket parameters
model_inputs.XCG_dry = rocket.COM_table(end, 2);
model_inputs.LREF = rocket.LREF;
model_inputs.SREF = rocket.SREF;
model_inputs.m_table = rocket.m_table;
model_inputs.Tx = rocket.Tx;
model_inputs.COM_table = rocket.COM_table;
model_inputs.I = rocket.MOI_dry;
model_inputs.MOI_table = rocket.MOI_table;

model_inputs.Rocket = rocket; % transition above to this 

end

