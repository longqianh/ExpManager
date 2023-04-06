classdef ExpManager < handle
properties
    exp_date 
    exp_name
    exp_main_dir
    exp_save_dir
    exp_data_dir
    toolbox_dir
end

methods(Static)
    function mkdirs(filePath)
    % Make directory if the file path does not exist
    % Longqian Huang, 2022.3.2
    
        if ~exist(filePath,'dir')
            [supPath,~] = fileparts(filePath);
            if ~exist(supPath,'dir')
                mkdirs(supPath)
            end
            mkdir(filePath)
        end
    end

    function subdir=get_subdir(maindir,sublevel)
        if ispc
            subdir_split=regexp(maindir,'\','split');
        else
            subdir_split=regexp(maindir,'/','split');
        end
        subdir="";
        for i=1:length(subdir_split)-sublevel
            subdir=fullfile(subdir,subdir_split{i});
        end
    end
   
    function import_toolboxes(toolbox_dir,exp_toolbox)
%       Put toolboxes in the Project dir (with the same level of current project)
      
        for i=1:length(exp_toolbox)
            toolbox_path=fullfile(toolbox_dir,exp_toolbox(i));
            addpath(genpath(toolbox_path));
            fprintf("Add toolbox %s\n",toolbox_path);
        end
    end
end

methods
    function obj=ExpManager(exp_name,exp_toolbox)
        cur_proj_dir=ExpManager.get_subdir(pwd,1);
        sup_proj_dir=ExpManager.get_subdir(pwd,2);
        obj.toolbox_dir=sup_proj_dir;
        obj.import_toolbox("Exp-toolbox");

        obj.exp_name=exp_name;
        today_ = string(datetime('today','InputFormat','yyyy-MM-dd'));
        obj.exp_date = erase(today_,'-');
       
        obj.exp_main_dir=pwd;
        obj.exp_save_dir=fullfile(sup_proj_dir,'Experiments',obj.exp_date,obj.exp_name);
        FileUtil.mkdirs(obj.exp_save_dir);
        obj.exp_data_dir=fullfile(cur_proj_dir,'data');
        obj.toolbox_dir=sup_proj_dir;
        ExpManager.import_toolboxes(sup_proj_dir,exp_toolbox);
    end

    function info(obj)
        fprintf('Exp Main dir: %s\n',obj.exp_main_dir);
        fprintf('Exp Data dir: %s\n',obj.exp_data_dir);
        fprintf('Exp Save dir: %s\n',obj.exp_save_dir);
    end

    function import_toolbox(obj,toolbox)
        ExpManager.import_toolboxes(obj.toolbox_dir,string(toolbox));
    end

   % todo: exp save wrapper
   
end
end