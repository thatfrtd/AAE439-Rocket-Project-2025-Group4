function [output] = missiledatcomimport_batch(filename, case_number, NMACH, NALPHA)

%% Initialization 
%%Initialize coeff 
initialsize = zeros(case_number, NALPHA, NMACH);
output.CD = initialsize;
output.CL = initialsize; % Lift coeff
output.XCP = initialsize; % Center of pressure position in calibers, measured from the
% moment reference center, divided by reference length. Positive values 

% Body axis static derivatives 
output.CA = initialsize; % Axial force coefficient
output.CY = initialsize; % Side force coefficient
output.CN = initialsize; % Normal force coefficient
output.CLL = initialsize; % Rolling moment coefficient
output.CM = initialsize; % Pitching moment coefficient
output.CLN = initialsize; % Yawing moment coefficient

output.CLLB = initialsize; % Rolling moment coefficient derivative with BETA
output.CYB = initialsize; % Side force coefficient derivative with BETA
output.CLNB = initialsize; % Yawing moment coefficient derivative with BETA

output.CNA = initialsize; % Normal force coefficient derivative with ALPHA
output.CMA = initialsize; % Pitching moment coefficient derivative with ALPHA

% Stability Derivatives
output.CAQ =  initialsize; % Axial force coefficient due to pitch rate
output.CYR = initialsize; % Side force coefficient due to yaw rate
output.CYP = initialsize; % Side force coefficient due to roll rate
output.CNQ = initialsize; % Normal force coefficient due to pitch rate
output.CNAD = initialsize; % Normal force coefficient due to rate of change of angle of attack

output.CLLP = initialsize; % Rolling moment coefficient due to roll rate
output.CLLR = initialsize; % Rolling moment coefficient due to yaw rate
output.CMQ = initialsize; % Pitching moment coefficient due to pitch rate
output.CMAD = initialsize; % Pitching moment coefficient due to rate of change of angle of attack
output.CLNP = initialsize; % Yawing moment coefficient due to roll rate
output.CLNR = initialsize; % Yawing moment coefficient due to yaw rate

%% ________________________
%% Open For006 File to read
fid2 = fopen(filename,"r");

%% _________________
%% Read Coefficients 
for c = 1:case_number
    %fprintf("%d\n",c)

    %Skip to Correct case
    case_string = string(c);
    if c < 10
        case_string = "  " + case_string;
    elseif c < 100
        case_string = " " + case_string;
    end
    indicator = 0;
    while indicator == 0
        line = fgetl(fid2);
        indicator = strcmp(line, sprintf("1         ***** THE USAF AUTOMATED MISSILE DATCOM * REV 3/99 *****     CASE %s", case_string));
    end

    ibool = 0;
    %Skip Header/Input Lines
    while ibool == 0
        line = fgetl(fid2);
        indicator = regexp(line,' {9}ALPHA','match');
        ibool = length(indicator);
    end

    for index = 1:NMACH
        StaticCoeffimport1 = fscanf(fid2,"%f");
     
        nrows = length(StaticCoeffimport1)/7;
    
        StaticCoeffimport1 = reshape(StaticCoeffimport1, [7 nrows]);
        StaticCoeff1 = transpose(StaticCoeffimport1);
        
        fgetl(fid2); %Skip Header Line
           
        
        StaticCoeffimport2 = fscanf(fid2,'%f');
        nrows2 = length(StaticCoeffimport2)/5;
        try
            StaticCoeffimport2 = reshape(StaticCoeffimport2, [5 nrows2]);   
        catch   %catch case does not account for another error case later in the dataset!!
            failline = fgetl(fid2);
            failure = regexp(failline,'\*+','match');
            failreplace = zeros(length(failure),1);
            remainder = regexp(failline, '\ +\d+\.\d+| +\-\d+\.\d+', 'match');
            remainderfloats = transpose(str2double(remainder));
            StaticCoeffimport2 = [StaticCoeffimport2; failreplace; remainderfloats];
            nrows2 = length(StaticCoeffimport2)/5;
            rowsleft = NALPHA - nrows2;
            if rowsleft > 0 
                remainder2 = [];
                for r = 1:rowsleft
                    dataline = fgetl(fid2);
                    remainderstr = regexp(dataline,'\ +\d+\.\d+| +\-\d+\.\d+','match');
                    remainder2floats = transpose(str2double(remainderstr));
                    remainder2 = [remainder2; remainder2floats];
                end
                StaticCoeffimport2 = [StaticCoeffimport2; remainder2];
                nrows2 = length(StaticCoeffimport2)/5;
            end
            
            StaticCoeffimport2 = reshape(StaticCoeffimport2, [5 nrows2]);  
        end   
        StaticCoeff2 = transpose(StaticCoeffimport2);
        
        ibool = 0;
        while ibool == 0
            line = fgetl(fid2);
            indicator = regexp(line,' {9}ALPHA','match');
            ibool = length(indicator);
        end
            
        DynamicCoeffimport1 = fscanf(fid2,'%f');
        nrows3 = length(DynamicCoeffimport1)/6;
        DynamicCoeffimport1 = reshape(DynamicCoeffimport1, [6 nrows3]);
        DynamicCoeff1 = transpose(DynamicCoeffimport1);
        
        ibool = 0;
        while ibool == 0
            line = fgetl(fid2);
            indicator = regexp(line,' {9}ALPHA','match');
            ibool = length(indicator);
        end
        
        DynamicCoeffimport2 = fscanf(fid2,'%f');
        nrows4 = length(DynamicCoeffimport2)/6;
        DynamicCoeffimport2 = reshape(DynamicCoeffimport2, [6 nrows4]);
        DynamicCoeff2 = transpose(DynamicCoeffimport2);
        
        ibool = 0;
        while ibool == 0
            line = fgetl(fid2);
            indicator = regexp(line,' {9}ALPHA','match');
            ibool = length(indicator);
        end

        DynamicCoeffimport3 = fscanf(fid2, '%f') ;
        nrows5 = length(DynamicCoeffimport3)/7;
        DynamicCoeffimport3 = reshape(DynamicCoeffimport3, [7 nrows5]);
        DynamicCoeff3 = transpose(DynamicCoeffimport3);

    %% ___________________
    %% Assign Coefficients 
        output.CN(c,:,index) = StaticCoeff1(:,2);
        output.CM(c,:,index) = StaticCoeff1(:,3);
        output.CA(c,:,index) = StaticCoeff1(:,4);
        output.CY(c,:,index) = StaticCoeff1(:,5);
        output.CLN(c,:,index) = StaticCoeff1(:,6);
        output.CLL(c,:,index) = StaticCoeff1(:,7);
    
        output.CL(c,:,index) = StaticCoeff2(:,2);
        output.CD(c,:,index) = StaticCoeff2(:,3);
    %     output.CL_CD(c,:,index) = StaticCoeff2(:,4);
        output.XCP(c,:,index) = StaticCoeff2(:,5);
        
        output.CNA(c,:,index) = DynamicCoeff1(:,2);
        output.CMA(c,:,index) = DynamicCoeff1(:,3);
        output.CYB(c,:,index) = DynamicCoeff1(:,4);
        output.CLNB(c,:,index) = DynamicCoeff1(:,5);
        output.CLLB(c,:,index) = DynamicCoeff1(:,6);

        output.CNQ(c,:,index) = DynamicCoeff2(:,2);
        output.CMQ(c,:,index) = DynamicCoeff2(:,3);
        output.CAQ(c,:,index) = DynamicCoeff2(:,4);
        output.CNAD(c,:,index) = DynamicCoeff2(:,5);
        output.CMAD(c,:,index) = DynamicCoeff2(:,6);
    
        output.CYR(c,:,index) = DynamicCoeff3(:,2);
        output.CLNR(c,:,index) = DynamicCoeff3(:,3);
        output.CLLR(c,:,index) = DynamicCoeff3(:,4);
        output.CYP(c,:,index) = DynamicCoeff3(:,5);
        output.CLNP(c,:,index) = DynamicCoeff3(:,6);
        output.CLLP(c,:,index) = DynamicCoeff3(:,7);

    %% __________________________
    %% Set Position for Next Loop
        if(index < NMACH)
            ibool = 0;
            while ibool == 0
                line = fgetl(fid2);
                indicator = regexp(line,' {9}ALPHA','match');
                ibool = length(indicator);
            end
        end
    end
end
end