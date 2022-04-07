mouse = 'FgDA_01';
protocol = 'Selina_C5D5R3E5R3';
filedir = ['D:/PhD/Photometry/DATA/round_20220307/' mouse '/' protocol '/Session Data/'] ;
DoricStudioVersion = '5.4.1.23';
cd(filedir)
files=uigetfile('*.mat','Select the INPUT DATA FILE(s) from BPOD','MultiSelect','on');
for i = 1:length(files)
    filename = files(i);
    sync_bpod_doric_data(filedir, filename,DoricStudioVersion)
end