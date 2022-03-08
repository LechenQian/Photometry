function MCoorTBL = table_MouseImplantFibers_MSFP

MCoorTBL = table('Size',[0 9],'VariableTypes',{'double','double','double',...
    'double','double','double','string','string','string'},'VariableNames',...
    {'Pos_HD19','Pos_Custom','Pos_Softw','AP_target','ML_target','DV_target'...
    'BrainRegion', 'ImplantHemisphere','VirusTransg'});

%                Standard       Custom
  % Positions in: HD19   B283-2001 HDFC-0002  DoricStudionROIindex  targetAP   targetML   targetDV    BrainRegion   Hemisphere side with implant   Virus or Transgenic
MCoorTBL(1,:) =  { 1,             1,                   0,             1.65,     0.60,       -4.0,     'NAccMedSh',   'L',                            'VT'};
MCoorTBL(2,:) =  { 3,             2,                   1,             1.25,     0.90,       -4.2,     'NAccCm',      'L',                            'VT'};
MCoorTBL(3,:) =  { 2,             3,                   2,             1.60,     1.25,       -4.5,     'NAccC',       'L',                            'VT'};
MCoorTBL(4,:) =  { 4,             4,                   3,             1.20,     1.75,       -4.5,     'NAccLatSh',   'L',                            'VT'};
MCoorTBL(5,:) =  { 5,             5,                   4,             0.70,     1.70,       -2.3,     'DMSa',        'L',                            'VT'};
MCoorTBL(6,:) =  { 6,             6,                   5,             0.20,     2.00,       -2.3,     'DMSm',        'L',                            'VT'};
MCoorTBL(7,:) =  { 7,             7,                   6,            -0.30,     2.30,       -2.3,     'DMSp',        'L',                            'VT'};
MCoorTBL(8,:) =  { 8,             8,                   7,             1.00,     2.30,       -2.3,     'DLSa',        'L',                            'VT'};
MCoorTBL(9,:) =  { 9,             9,                   8,             0.45,     2.70,       -2.3,     'DLSm',        'L',                            'VT'};
MCoorTBL(10,:) = { 10,            10,                  9,             0.00,     3.00,       -2.3,     'DLSp',        'L',                            'VT'}; 
MCoorTBL(11,:) = { 15,            11,                  10,           -0.80,     2.75,       -3.0,     'TSam',        'L',                            'VT'};
MCoorTBL(12,:) = { 16,            12,                  11,           -1.30,     3.00,       -2.0,     'TSal',        'L',                            'VT'};
MCoorTBL(13,:) = { 18,            13,                  12,           -0.60,     3.40,       -3.0,     'TS',          'L',                            'VT'}; 
MCoorTBL(14,:) = { 19,            14,                  13,           -1.40,     3.65,       -3.3,     'TSp',         'L',                            'VT'};
MCoorTBL(15,:) = { nan,           15,                  14,             nan,      nan,        nan,     'ctr',         'L',                            'VT'}; 

brainloc = MCoorTBL{:,'BrainRegion'};
MCoorTBL.Properties.RowNames = cellstr(brainloc); 

% To access coordinates of a mouse:
% MCoorTBL('SM106','AP_Coord_center') 
% table2array(MCoorTBL('SM106',[3 5 6]))

end