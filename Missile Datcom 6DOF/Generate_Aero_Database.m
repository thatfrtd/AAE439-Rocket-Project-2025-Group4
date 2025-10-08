function [aerotable] = Generate_Aero_Database(ALT, PHI, ALPHA, MACH, SREF, LREF, XCG_dry, NC_type, NC_power, NC_length,NC_diam, AF1_length,AF1_diam, FN1_number, FN1_type, FN1_NACA, FN1_height, FN1_root, FN1_tip, FN1_XLE, FN1_sweep, FN1_edgeRad, ZUPPER, ZLOWER, LMAXU, LMAXL, LFLATU, LFLATL, ROUGH, DEXIT, target_batch_size)
%GENERATE_AERO_DATABASE Loop through flight conditions and call DATCOM to create aerotable
%   For each ALT and PHI, calls MissleDat_Sustainer_opt which writes to the 
% input files and runs missledatcom.exe. 

cd("Missile Datcom 6DOF")

%% 4D Table Generation
NALT = size(ALT, 2); % Number of Altitudes 
NPHI = size(PHI, 2); % Number of Total Aerodynamic Clock Angles
NALPHA = size(ALPHA, 2); % Number of AoA 
NMACH = size(MACH, 2); % Number of Mach numbers

database_initial = zeros(NALT,NPHI,NALPHA,NMACH);

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

% Initialize loading bar
loading_bar = waitbar(0, 'Aerotable Progress: 0 %');
iteration = 0;
n = NALT * NPHI;

% Calculate number of batches (and DATCOM calls) needed
num_batch = (ceil(n/target_batch_size));

% Altitude loop
for b = 1:num_batch
    % Find the cases that are in this batch
    if b == num_batch
        batch_cases = ((b - 1) * target_batch_size + 1) : n;
    else
        batch_cases = (((b - 1) * target_batch_size + 1) : (b * target_batch_size));
    end

    % Calculate flight alt and phi for each case in batch
    n_alt = ceil(batch_cases / NPHI);
    n_phi = mod(batch_cases - 1, NPHI) + 1;

    [CA,CY,CN,CLL,CM,CLN,CMQ,CMAD,CLP,CLR,CNP,CNR,CAQ,CNQ,CNAD,CYR,CYP,CLB,CYB,CNB,XCP,CD,CL] = ...
            Run_Missile_Datcom(SREF, LREF, XCG_dry, ALT(n_alt), PHI(n_phi), NMACH, MACH, ...
                NALPHA, ALPHA,NC_type,NC_power,NC_length,NC_diam,AF1_length,AF1_diam, ...
                FN1_number, FN1_type, FN1_NACA, FN1_height, ...
                FN1_root, FN1_tip, FN1_XLE, FN1_sweep, FN1_edgeRad, ...
                ZUPPER, ZLOWER, LMAXU,LMAXL,LFLATU,LFLATL, ROUGH, DEXIT, "");

    for c = 1:length(batch_cases)
        i = n_alt(c);
        j = n_phi(c);
        
        aerotable.CA_sustainer(i,j,:,:) = CA(c,:,:);
        aerotable.CY_sustainer(i,j,:,:) = CY(c,:,:);
        aerotable.CN_sustainer(i,j,:,:) = CN(c,:,:);
        aerotable.CLL_sustainer(i,j,:,:) = CLL(c,:,:);
        aerotable.CM_sustainer(i,j,:,:) = CM(c,:,:); 
        aerotable.CLN_sustainer(i,j,:,:) = CLN(c,:,:);
        aerotable.CMQ_sustainer(i,j,:,:) = CMQ(c,:,:);
        aerotable.CMAD_sustainer(i,j,:,:) = CMAD(c,:,:);
        aerotable.CLP_sustainer(i,j,:,:) = CLP(c,:,:);
        aerotable.CLR_sustainer(i,j,:,:) = CLR(c,:,:);
        aerotable.CNP_sustainer(i,j,:,:) = CNP(c,:,:);
        aerotable.CNR_sustainer(i,j,:,:) = CNR(c,:,:);
        aerotable.CAQ_sustainer(i,j,:,:) = CAQ(c,:,:);
        aerotable.CNQ_sustainer(i,j,:,:) = CNQ(c,:,:);
        aerotable.CNAD_sustainer(i,j,:,:) = CNAD(c,:,:);
        aerotable.CYR_sustainer(i,j,:,:) = CYR(c,:,:);
        aerotable.CYP_sustainer(i,j,:,:) = CYP(c,:,:);
        aerotable.CLB_sustainer(i,j,:,:) = CLB(c,:,:);
        aerotable.CYB_sustainer(i,j,:,:) = CYB(c,:,:); 
        aerotable.CNB_sustainer(i,j,:,:) = CNB(c,:,:);
        aerotable.XCP_sustainer(i,j,:,:) = XCP(c,:,:);
        aerotable.CD_sustainer(i,j,:,:) = CD(c,:,:);
        aerotable.CL_sustainer(i,j,:,:) = CL(c,:,:);
    end

    % Update loading bar
    iteration = iteration + length(batch_cases);
    waitbar(iteration/n, loading_bar, sprintf('Aerotable Progress: %d %%', floor(iteration/n*100)));
end

% Close loading bar
close(loading_bar)

cd("../")

end

