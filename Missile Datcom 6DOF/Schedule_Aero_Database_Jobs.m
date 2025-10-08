function [aerotable] = Schedule_Aero_Database_Jobs(ALT, PHI, ALPHA, MACH, SREF, LREF, XCG_dry, NC_type, NC_power, NC_length,NC_diam, AF1_length,AF1_diam, FN1_number, FN1_type, FN1_NACA, FN1_height, FN1_root, FN1_tip, FN1_XLE, FN1_sweep, FN1_edgeRad, ZUPPER, ZLOWER, LMAXU, LMAXL, LFLATU, LFLATL, ROUGH, DEXIT, target_batch_size, thread_number)
%SCHEDULE_AERO_DATABASE_JOBS Schedule batch jobs accross multiple threads to generate aero database 
%   Schedules multiple Generate_Aero_Database_Parallel to run on seperate
%   threads to multithread table generation.

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
loading_bar = waitbar(0, 'Creating Aerotable Jobs: 0 %');
iteration = 0;
n = NALT * NPHI;

% Calculate number of batches (and DATCOM calls) needed
thread_case_number = (ceil(n/thread_number));

% Create cluster object to run jobs on
if any(string(parallel.listProfiles) == "AerotableBatch") == false
    parallel.importProfile('Missile Datcom 6DOF/AerotableBatch.mlsettings');
end
clust = parcluster("AerotableBatch");

% Initialize jobs and job folders
for t = 1:thread_number
    % Create jobs so that tasks can be added to them
    jobs(t) = createJob(clust);

    % Create folder for job if needed
    if ~exist(string("BatchJobFolders/Folder" + t), 'dir')
        copyfile("BatchJobFolders/Folder1", string("BatchJobFolders/Folder" + t));
    end
end

for t = 1:thread_number
    % Find the cases that are in this thread
    if t == thread_number
        thread_cases = ((t - 1) * thread_case_number + 1) : n;
    else
        thread_cases = (((t - 1) * thread_case_number + 1) : (t * thread_case_number));
    end

    % Calculate flight alt and phi for each case in batch
    n_alt = ceil(thread_cases / NPHI);
    n_phi = mod(thread_cases - 1, NPHI) + 1;
    
    tasks(t) = createTask(jobs(t), @Generate_Aero_Database_Job, 1, {ALT(n_alt), PHI(n_phi), ALPHA, MACH, SREF, LREF, XCG_dry, NC_type, ...
        NC_power, NC_length,NC_diam, AF1_length,AF1_diam, FN1_number, FN1_type, FN1_NACA, FN1_height, ...
        FN1_root, FN1_tip, FN1_XLE, FN1_sweep, FN1_edgeRad, ZUPPER, ZLOWER, LMAXU, LMAXL, LFLATU, ...
        LFLATL ROUGH, DEXIT, target_batch_size, string("BatchJobFolders/Folder" + t)},'CaptureDiary',false);
    
    % Update loading bar
    iteration = iteration + length(thread_cases);
    waitbar(iteration/n, loading_bar, sprintf('Creating Aerotable Jobs: %d %%', floor(iteration/n*100)));
end

% Submit Jobs
for t = 1:thread_number
    submit(jobs(t))
end

% Look at diary of job for real-time updates? jobs(1).diary shows number of
% batches completed by job one at the time of calling
waitbar(0, loading_bar, sprintf('Aerotable Progress: %d %%', 0));
for t = 1:thread_number
    wait(jobs(t));

    waitbar(t/thread_number, loading_bar, sprintf('Aerotable Progress: %d %%', floor(t/thread_number*100)));
end

% Close loading bar
close(loading_bar)

for t = 1:thread_number
    % Find the cases that are in this thread
    if t == thread_number
        thread_cases = ((t - 1) * thread_case_number + 1) : n;
    else
        thread_cases = (((t - 1) * thread_case_number + 1) : (t * thread_case_number));
    end

    % Calculate flight alt and phi for each case in batch
    n_alt = ceil(thread_cases / NPHI);
    n_phi = mod(thread_cases - 1, NPHI) + 1;
    
    batch_table = fetchOutputs(jobs(t));

    for c = 1:length(thread_cases)
        i = n_alt(c);
        j = n_phi(c);

        aerotable.CA_sustainer(i,j,:,:) = batch_table{1}.CA_sustainer(c,:,:);
        aerotable.CY_sustainer(i,j,:,:) = batch_table{1}.CY_sustainer(c,:,:);
        aerotable.CN_sustainer(i,j,:,:) = batch_table{1}.CN_sustainer(c,:,:);
        aerotable.CLL_sustainer(i,j,:,:) = batch_table{1}.CLL_sustainer(c,:,:);
        aerotable.CM_sustainer(i,j,:,:) = batch_table{1}.CM_sustainer(c,:,:); 
        aerotable.CLN_sustainer(i,j,:,:) = batch_table{1}.CLN_sustainer(c,:,:);
        aerotable.CMQ_sustainer(i,j,:,:) = batch_table{1}.CMQ_sustainer(c,:,:);
        aerotable.CMAD_sustainer(i,j,:,:) = batch_table{1}.CMAD_sustainer(c,:,:);
        aerotable.CLP_sustainer(i,j,:,:) = batch_table{1}.CLP_sustainer(c,:,:);
        aerotable.CLR_sustainer(i,j,:,:) = batch_table{1}.CLR_sustainer(c,:,:);
        aerotable.CNP_sustainer(i,j,:,:) = batch_table{1}.CNP_sustainer(c,:,:);
        aerotable.CNR_sustainer(i,j,:,:) = batch_table{1}.CNR_sustainer(c,:,:);
        aerotable.CAQ_sustainer(i,j,:,:) = batch_table{1}.CAQ_sustainer(c,:,:);
        aerotable.CNQ_sustainer(i,j,:,:) = batch_table{1}.CNQ_sustainer(c,:,:);
        aerotable.CNAD_sustainer(i,j,:,:) = batch_table{1}.CNAD_sustainer(c,:,:);
        aerotable.CYR_sustainer(i,j,:,:) = batch_table{1}.CYR_sustainer(c,:,:);
        aerotable.CYP_sustainer(i,j,:,:) = batch_table{1}.CYP_sustainer(c,:,:);
        aerotable.CLB_sustainer(i,j,:,:) = batch_table{1}.CLB_sustainer(c,:,:);
        aerotable.CYB_sustainer(i,j,:,:) = batch_table{1}.CYB_sustainer(c,:,:); 
        aerotable.CNB_sustainer(i,j,:,:) = batch_table{1}.CNB_sustainer(c,:,:);
        aerotable.XCP_sustainer(i,j,:,:) = batch_table{1}.XCP_sustainer(c,:,:);
        aerotable.CD_sustainer(i,j,:,:) = batch_table{1}.CD_sustainer(c,:,:);
        aerotable.CL_sustainer(i,j,:,:) = batch_table{1}.CL_sustainer(c,:,:);
    end
end

delete(jobs)

end

