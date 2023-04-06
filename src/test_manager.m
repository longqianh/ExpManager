clc;close all;clear;

exp_toolbox=["SLM-toolbox","AO-toolbox"];
ma=ExpManager('exp_test',exp_toolbox);
ma.info()

%% single toolbox import
ma.import_toolbox('QPSL');

%% exp data saving
test_data=rand(10);
data_name='test';
save(fullfile(ma.exp_save_dir,data_name),"test_data");
