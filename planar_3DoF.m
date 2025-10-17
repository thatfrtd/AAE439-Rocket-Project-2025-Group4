function xdot = rail_sim_derivative(t, x, rocket_params, launch_site)
    % Unpack state
    r_x = x(1); r_y = x(2); % position in inertial frame
    v_x = x(3); v_y = x(4); % velocities in inertial frame
    psi = x(5); w = x(6);
    
    g = 9.81; % [m / s2]
    [~, a_sound, ~, rho] = atmoscoesa(r_y + launch_site.elevation); % [Temperature, speed o' sound, pressure, air density]

    mass =  interp1(rocket_params.mass_table(:, 1), rocket_params.mass_table(:, 2), t, "linear", "extrap");
    thrust = interp1(rocket_params.thrust_table(:, 1), rocket_params.thrust_table(:, 2), t, "linear", 0); % thrust from table + interpolation

    V_sqr = v_x ^ 2 + v_y ^ 2;
    Q = 1 / 2 * V_sqr * rho;

    a_friction = -cosd(launch_site.elevation_angle) * launch_site.rail_mu * g;
    a_gravity = -sind(launch_site.elevation_angle) * g;
    F_drag = -Q * rocket_params.axial_coeff(0) * rocket_params.cross_area;
    F_thrust = thrust;

    dvdt_mag = a_friction + a_gravity + (F_drag + F_thrust) / mass; 
    dvdt_direction = [cosd(launch_site.elevation_angle); sind(launch_site.elevation_angle)];

    dvdt = dvdt_mag * dvdt_direction;

    xdot = [v_x; v_y; dvdt; 0; 0];
end

function [raily_left, isterminal, direction] = rail_event_function(t, x, rocket_params, launch_site)
    raily_left = sind(launch_site.elevation_angle) * launch_site.rail_length - x(2);
    isterminal = 1; % Stop simulation
    direction = -1; % Trigger when zero is approached from positive side
end

function dstatedt = planar_3dof(t, state, rocket_params, launch_site, wind)
    % Unpack state
    r_x = state(1); r_y = state(2); % position in inertial frame
    v_x = state(3); v_y = state(4); % velocities in inertial frame
    psi = state(5); w = state(6);
    mass =  interp1(rocket_params.mass_table(:, 1), rocket_params.mass_table(:, 2), t, "linear", "extrap");
    thrust = interp1(rocket_params.thrust_table(:, 1), rocket_params.thrust_table(:, 2), t, "linear", 0); % thrust from table + interpolation
    COM = interp1(rocket_params.COM_table(:, 1), rocket_params.COM_table(:, 2), t, "linear", "extrap");
    % Unpack parameters
    Iz = rocket_params.inertia * mass; % Moment of Inertia of the body
    g = 9.81; % [m / s2]
    gvec = [0; -g];
    R = [cos(psi), -sin(psi); sin(psi), cos(psi)]; % DCM from body to inertial frame
    velocity_b = R' * [v_x; v_y];
    alpha = atan2(velocity_b(2), velocity_b(1));
    % Drag claculation based on equation
    [~, a_sound, ~, rho] = atmoscoesa(r_y + launch_site.elevation); % [Temperature, speed o' sound, pressure, air density]
    Q = 0.5 * rho * norm(velocity_b)^2;
    A_axial = Q * rocket_params.cross_area * rocket_params.axial_coeff(alpha);
    A_normal = Q * rocket_params.cross_area * rocket_params.normal_coeff(alpha);
    F_Drag = [-A_axial;
              -A_normal];
    Moment_aero = Q * rocket_params.b_ref * rocket_params.cross_area * rocket_params.pitch_moment(alpha);
    Moment_aeroforce = cross([(COM - rocket_params.COP(alpha)); 0; 0],  [F_Drag; 0]);
    Moment = Moment_aero + Moment_aeroforce(3);
    F_Thrust = [thrust; 0];
    F_Net = F_Drag + F_Thrust; % body frame forces
    % Derivative calculations
    posdot = [v_x; v_y]; % velocity in inertial frame
    psidot = w; % angular velocity relation
    dvdt = R*(F_Net/mass) + gvec; % - R*[-w*velocity_b(2); w*velocity_b(1)];
    wdot = Moment/Iz; % torque = I*alpha
    % mass change rate based on Thrust
    dstatedt = [posdot; dvdt; psidot; wdot];
end

function [vy, isterminal, direction] = apogee_event_function(t, x)
    vy = x(4);
    isterminal = 1; % Stop simulation
    direction = -1; % Trigger when zero is approached from positive side
end

rocket_file_name = "hitech";
rocket = load("Rockets\" + rocket_file_name + ".mat").rocket;
% Parameters
rocket_params.radius = rocket.AF1_diam/2; % [m]
rocket_params.length = rocket.AF1_length + rocket.NC_length; % [m]
rocket_params.mass_table = rocket.motor.mass_table;  % kg
rocket_params.mass_table = [[0, rocket_params.mass_table(1, 2)]; rocket_params.mass_table; [rocket_params.mass_table(end, 1) + 1e-3, rocket_params.mass_table(end, 2)]];
rocket_params.thrust_table = rocket.motor.thrust_table;
rocket_params.COM_table = rocket.COM_table; % [m] from nosecone
rocket_params.COM_table = [[0, rocket_params.COM_table(1, 2)]; rocket_params.COM_table; [rocket_params.COM_table(end, 1) + 1e-3, rocket_params.COM_table(end, 2)]];
rocket_params.cross_area = rocket.SREF; % m^2
rocket_params.axial_coeff = @(alpha) 0.475;
rocket_params.normal_coeff = @(alpha) 0 + 2 * alpha;
rocket_params.pitch_moment = @(alpha) alpha * 2.75;
rocket_params.b_ref = rocket.AF1_diam; % [m]
rocket_params.inertia = (rocket_params.radius ^ 2 / 4 + rocket_params.length ^ 2 / 12); % huh?
rocket_params.Isp = 200; % [s]
rocket_params.COP = @(alpha) 1.02;%(1 - 0.1 * (rad2deg(abs(alpha)) / 10)) * 1.5; % [m] from nosecone

launch_site = load("LaunchSites\ThunderStruck.mat").launch_site;
launch_site.rail_length = 1.5;

wind = load("Winds\ConstantNoWind.mat").wind;
wind.wind_speed = 4; % [m / s]

%%
% Initial state: [x, y, u, v, psi, w]
x0_rail = [0; 0; 0; 0; deg2rad(launch_site.elevation_angle); deg2rad(0)];
% Time span
tspan = [0 9];

%% Rail Sim
% Set event functions for stopping conditions, etc
odeoptions = odeset(Events = @(t, x) rail_event_function(t, x, rocket_params, launch_site));

% Run ODE
[t_rail, x_rail] = ode45(@(t, state) rail_sim_derivative(t, state, rocket_params, launch_site), tspan, x0_rail, odeoptions);

x0_3DoF = x_rail(end, :);

%% Set event functions for stopping conditions, etc
odeoptions = odeset(Events = @apogee_event_function);

% Run ODE
[t, state] = ode45(@(t, state) planar_3dof(t, state, rocket_params, launch_site, wind), [t_rail(end), tspan(end)], x0_3DoF, odeoptions);
% ode only puts in t, state. We need @(t,state) to act as a adapater for
% ode. Then we can add in params with ode "thinking" it only put in
% (t,state)
% The solver is the one that is determining the states based on the
% previous state/initial data.
% The translational_3dof is handling getting the velocity and acceleration
% of the current state (F=ma and grabbing u, v, and w).
% ode45 uses an estimation method with error checking to get as closec as
% possible to the correct next state
% Terminal Velocity
% Calculate terminal velocity
state = state';
velocity_b = zeros(2, numel(t));
alpha = zeros(1, numel(t)); % angle of attack
gamma = zeros(1, numel(t)); % flightpath angle
for i = 1:numel(t)
    velocity = state(3:4, i);
    psi = state(5, i);
    R = [cos(psi), sin(psi); -sin(psi), cos(psi)]; % DCM inertial to body
    velocity_b(1:2, i) = R * velocity; % convert from inertial to body
    alpha(i) = atan2(velocity_b(2,i), velocity_b(1,i)); % angle of attack calculation
    gamma(i) = atan2(-velocity(2), velocity(1));
end
figure
tiledlayout(2, 2)
nexttile
plot(t, state(1, :)); hold on
plot(t, state(2, :)); hold off
xlabel("Time [s]")
ylabel("Position [m]")
legend("x", "y")
nexttile
plot(t, state(3, :)); hold on
plot(t, state(4, :)); hold off
xlabel("Time [s]")
ylabel("Velocity [m / s]")
legend("u", "v")
nexttile
plot(t, wrapTo180(rad2deg(state(5, :))));
xlabel("Time [s]")
ylabel("Orientation [deg]")
nexttile
plot(t, rad2deg(state(6, :)));
xlabel("Time [s]")
ylabel("Angular Velocity [deg / s]")
figure
plot(state(1, :), state(2,:), "-");
title("Trajectory")
xlabel("X [m]")
ylabel("Y [m]")
axis equal
grid on
figure
plot(t, rad2deg(alpha));
xlabel("Time [s]")
ylabel("Angle of Attack [deg]")
figure
plot(t, rad2deg(gamma));
xlabel("Time [s]")
ylabel("Flight Path Angle [deg]")
% Plot