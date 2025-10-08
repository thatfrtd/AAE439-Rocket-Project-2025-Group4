classdef ParachuteDeploymentType < uint32
    %PARACHUTEDEPLOYMENTTYPE Summary of this class goes here
    %   Detailed explanation goes here
    %   MUST BE FILE NAMES OF PARACHUTE DEPLOYMENT MODEL SUBSYSTEM FILES
    enumeration
        Apogee_Delay_Deployment (1)
        Burnout_Delay_Deployment (2)
        Elevation_Deployment (3)
        Custom_Deployment_C (4)
        Custom_Deployment_Matlab (5)
        Trigger_Signal_Deployment (6)
    end
end


