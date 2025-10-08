function [aerotable] = Generate_Aero_Database_Job(ALT, PHI, ALPHA, MACH, SREF, LREF, XCG_dry, NC_type, NC_power, NC_length,NC_diam, AF1_length,AF1_diam, FN1_number, FN1_type, FN1_NACA, FN1_height, FN1_root, FN1_tip, FN1_XLE, FN1_sweep, FN1_edgeRad, ZUPPER, ZLOWER, LMAXU, LMAXL, LFLATU, LFLATL, ROUGH, DEXIT, target_batch_size, CurrentFolder)
%GENERATE_AERO_DATABASE_JOB Loop through flight conditions and call DATCOM to create aerotable
%   For each ALT and PHI, calls Run_Missile_Datcom which writes to the 
% input files and runs missledatcom.exe. Specifically tailored to be called
% by Schedule_Aero_Database_Jobs and run on seperate threads.

cd(CurrentFolder);

%% 4D Table Generation
case_number = size(ALT, 2);
NALPHA = size(ALPHA, 2); % Number of AoA 
NMACH = size(MACH, 2); % Number of Mach numbers

database_initial = zeros(case_number,NALPHA,NMACH);

% Sustainer Initialization
aerotable = struct();
aerotable.CA_sustainer = database_initial;
aerotable.CY_sustainer = database_initial;
aerotable.CN_sustainer = database_initial;
aerotable.CLL_sustainer = database_initial;
aerotable.CM_sustainer = database_initial; 
aerotable.CLN_sustainer = database_initial; 
aerotable.CMQ_sustainer = database_initial;
aerotable.CMAD_sustainer = database_initial;
aerotable.CLP_sustainer = database_initial;
aerotable.CLR_sustainer = database_initial;
aerotable.CNP_sustainer = database_initial; 
aerotable.CNR_sustainer = database_initial;
aerotable.CAQ_sustainer = database_initial;
aerotable.CNQ_sustainer = database_initial;
aerotable.CNAD_sustainer = database_initial;
aerotable.CYR_sustainer = database_initial;
aerotable.CYP_sustainer = database_initial;
aerotable.CLB_sustainer = database_initial;
aerotable.CYB_sustainer = database_initial; 
aerotable.CNB_sustainer = database_initial;
aerotable.XCP_sustainer = database_initial;
aerotable.CD_sustainer = database_initial;
aerotable.CL_sustainer = database_initial; 

% Calculate number of batches (and DATCOM calls) needed
num_batch = (ceil(case_number/target_batch_size));

% Altitude loop
for b = 1:num_batch
    % Find the cases that are in this batch
    if b == num_batch
        batch_cases = ((b - 1) * target_batch_size + 1) : case_number;
    else
        batch_cases = (((b - 1) * target_batch_size + 1) : (b * target_batch_size));
    end

    [CA,CY,CN,CLL,CM,CLN,CMQ,CMAD,CLP,CLR,CNP,CNR,CAQ,CNQ,CNAD,CYR,CYP,CLB,CYB,CNB,XCP,CD,CL] = ...
            Run_Missile_Datcom(SREF, LREF, XCG_dry, ALT(batch_cases), PHI(batch_cases), NMACH, MACH, ...
                NALPHA, ALPHA,NC_type,NC_power,NC_length,NC_diam,AF1_length,AF1_diam, ...
                FN1_number, FN1_type, FN1_NACA, FN1_height, ...
                FN1_root, FN1_tip, FN1_XLE, FN1_sweep, FN1_edgeRad, ...
                ZUPPER, ZLOWER, LMAXU,LMAXL,LFLATU,LFLATL, ROUGH, DEXIT, "");

    aerotable.CA_sustainer(batch_cases, :, :) = CA;
    aerotable.CY_sustainer(batch_cases, :, :) = CY;
    aerotable.CN_sustainer(batch_cases, :, :) = CN;
    aerotable.CLL_sustainer(batch_cases, :, :) = CLL;
    aerotable.CM_sustainer(batch_cases, :, :) = CM; 
    aerotable.CLN_sustainer(batch_cases, :, :) = CLN;
    aerotable.CMQ_sustainer(batch_cases, :, :) = CMQ;
    aerotable.CMAD_sustainer(batch_cases, :, :) = CMAD;
    aerotable.CLP_sustainer(batch_cases, :, :) = CLP;
    aerotable.CLR_sustainer(batch_cases, :, :) = CLR;
    aerotable.CNP_sustainer(batch_cases, :, :) = CNP;
    aerotable.CNR_sustainer(batch_cases, :, :) = CNR;
    aerotable.CAQ_sustainer(batch_cases, :, :) = CAQ;
    aerotable.CNQ_sustainer(batch_cases, :, :) = CNQ;
    aerotable.CNAD_sustainer(batch_cases, :, :) = CNAD;
    aerotable.CYR_sustainer(batch_cases, :, :) = CYR;
    aerotable.CYP_sustainer(batch_cases, :, :) = CYP;
    aerotable.CLB_sustainer(batch_cases, :, :) = CLB;
    aerotable.CYB_sustainer(batch_cases, :, :) = CYB; 
    aerotable.CNB_sustainer(batch_cases, :, :) = CNB;
    aerotable.XCP_sustainer(batch_cases, :, :) = XCP;
    aerotable.CD_sustainer(batch_cases, :, :) = CD;
    aerotable.CL_sustainer(batch_cases, :, :) = CL;

    fprintf("Completed Batch %g\n", b)
end

end

