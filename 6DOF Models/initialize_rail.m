function [rail_inputs] = initialize_rail(rocket, launch_site)
%INITIALIZE_RAIL Summary of this function goes here
%   Detailed explanation goes here
arguments
    rocket {mustBeA(rocket, "Rocket")}
    launch_site {mustBeA(launch_site, "LaunchSite")}
end

%% Rail Sim
rail_inputs = struct();
rail_inputs.Cd = 0.6;
rail_inputs.e = deg2rad(launch_site.elevation_angle);
rail_inputs.l = launch_site.rail_length; % [m] launch rail 
rail_inputs.mu = launch_site.rail_mu;
rail_inputs.S = rocket.SREF * 0.0254^2;
rail_inputs.T_data = rocket.Tx(:, 2);
rail_inputs.T_time = rocket.Tx(:, 1);
rail_inputs.m_data = rocket.m_table(:, 2);
rail_inputs.m_time = rocket.m_table(:, 1);

end

