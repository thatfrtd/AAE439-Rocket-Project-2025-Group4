function [initial_conditions] = SixDoF_output_to_6DoF_initial_conditions(sixDoF_outputs, launch_site)
%SIXDOF_OUTPUT_TO_6DOF_INITIAL_CONDITIONS Summary of this function goes here
%   Detailed explanation goes here

initial_conditions = struct();
initial_conditions.time0 = sixDoF_outputs.time(end);
initial_conditions.HG0 = sixDoF_outputs.Altitude(end);
initial_conditions.x0_NED = sixDoF_outputs.x_NED(:, end);
initial_conditions.launch_elevation = launch_site.elevation;
initial_conditions.uvw0_B = sixDoF_outputs.uvw_B(:, end);
final_rotation = sixDoF_outputs.Rotation_Angles(end, :);
initial_conditions.BG0 = angle2dcm(final_rotation(1), final_rotation(2), final_rotation(3));
initial_conditions.lambda0 = sixDoF_outputs.lat(end);
initial_conditions.mu0 = sixDoF_outputs.long(end);
initial_conditions.w0_BG = sixDoF_outputs.w_BG(:, end);

end

