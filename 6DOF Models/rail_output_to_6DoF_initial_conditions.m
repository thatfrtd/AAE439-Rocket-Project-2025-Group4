function [initial_conditions] = rail_output_to_6DoF_initial_conditions(rail_outputs, launch_site)
%RAIL_OUTPUT_TO_6DOF_INITIAL_CONDITIONS Summary of this function goes here
%   Detailed explanation goes here
arguments
    rail_outputs
    launch_site {mustBeA(launch_site, "LaunchSite")}
end

initial_conditions = struct();

initial_conditions.time0 = rail_outputs.Time(end);

initial_conditions.HG0 = launch_site.elevation + rail_outputs.Height(end); % [m] Initial Altitude
initial_conditions.launch_elevation = launch_site.elevation;

initial_conditions.x0_NED = [0,0,-initial_conditions.HG0]; % [m] Initial Altitude

initial_conditions.uvw0_B = [rail_outputs.Velocity(end);0;0]; % [m / s]

rotAngles0 = [launch_site.compass_angle, ...
              launch_site.elevation_angle,0]; % Yaw, Pitch, Roll [deg]
rotAngles0 = deg2rad(rotAngles0); % Convert to radian for simulink
psi0 = rotAngles0(1);
theta0 = rotAngles0(2);
phi0 = rotAngles0(3);

% Body to Geodetic DCM
R_3 = [cos(psi0),sin(psi0),0;
       -sin(psi0),cos(psi0),0;
       0,0,1];
R_2 = [cos(theta0),0,-sin(theta0);
       0,1,0;
       sin(theta0),0,cos(theta0)];
R_1 = [1,0,0;
       0,cos(phi0),sin(phi0);
       0,-sin(phi0),cos(phi0)];

initial_conditions.BG0 = R_1*R_2*R_3; % Geodetic to Body

latlong0 = deg2rad([launch_site.latitude; ...
                    launch_site.longitude]);
initial_conditions.lambda0 = latlong0(1);
initial_conditions.mu0 = latlong0(2);

% Initial Spin Rate (Body Frame)
w0_BG = [0;0;0]; % [deg/sec] PQR (spin relative to geodetic)
initial_conditions.w0_BG = deg2rad(w0_BG);

end

