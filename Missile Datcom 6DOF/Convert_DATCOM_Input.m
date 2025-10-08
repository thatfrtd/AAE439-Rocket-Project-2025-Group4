function [SREF, LREF, XCG_dry, ALT, NC_length, NC_diam, AF1_length, AF1_diam, FN1_height, FN1_root, FN1_tip, FN1_XLE, FN1_edgeRad, ROUGH, DEXIT] = Convert_DATCOM_Input(SREF, LREF, XCG_dry, ALT, NC_length, NC_diam, AF1_length, AF1_diam, FN1_height, FN1_root, FN1_tip, FN1_XLE, FN1_edgeRad, ROUGH, DEXIT)
% Function converts variables from metric to imperial
%   Missile DATCOM needs inputs in imperial form. This function converts
%   all metric units into imperial for input into Missile DATCOM
SREF = 1550*SREF; % [m2] to [in2]
LREF = 39.3701*LREF; % [m] to [in]
XCG_dry = 39.3701*XCG_dry; % [m] to [in]
%ALT = 39.3701*ALT; % [m] to [in]
NC_length = 39.3701*NC_length; % [m] to [in]
NC_diam = 39.3701*NC_diam; % [m] to [in]
AF1_length = 39.3701*AF1_length; % [m] to [in]
AF1_diam = 39.3701*AF1_diam; % [m] to [in]
FN1_height = 39.3701*FN1_height; % [m] to [in]
FN1_root = 39.3701*FN1_root; % [m] to [in]
FN1_tip = 39.3701*FN1_tip; % [m] to [in]
FN1_XLE = 39.3701*FN1_XLE; % [m] to [in]
FN1_edgeRad = 39.3701*FN1_edgeRad; % [m] to [in]
ROUGH = 39.3701*ROUGH; % [m] to [in]
DEXIT = 39.3701*DEXIT; % [m] to [in]
end