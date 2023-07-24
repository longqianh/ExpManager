classdef ExpManager < handle
properties
    exp_date 
    exp_name
    exp_dir
    exp_main_dir
    exp_save_dir
    exp_data_dir
    toolbox_dir
end

methods(Static)
    function filePath=mkdirs(filePath)
    % Make directory if the file path does not exist
    % Longqian Huang, 2022.3.2
    
        if ~exist(filePath,'dir')
            [supPath,~] = fileparts(filePath);
            if ~exist(supPath,'dir')
                mkdir(supPath);
                fprintf("make sup dirs: %s\n", supPath);
            end
            mkdir(filePath);
            fprintf("make dir: %s\n", filePath);
        else
            fprintf("path already exists.\n");
        end
    end

    function subdir=get_subdir(maindir,sublevel)
        if ispc
            subdir_split=regexp(maindir,'\','split');
            subdir="";
        else
            subdir_split=regexp(maindir,'/','split');
            subdir="/";
        end
        
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
    function obj=ExpManager(exp_name,exp_toolbox,exp_date)
        % When exp_name == [], ExpManager will not mkdirs
        
        cur_proj_dir=ExpManager.get_subdir(pwd,1);
        sup_proj_dir=ExpManager.get_subdir(pwd,2);
        obj.toolbox_dir=sup_proj_dir;
        addpath(genpath(cur_proj_dir));
        obj.import_toolbox("Exp-toolbox");

        obj.exp_name=exp_name;
        if nargin<3
            today_ = string(datetime('today','InputFormat','yyyy-MM-dd'));
            obj.exp_date = erase(today_,'-');
        else
            obj.exp_date=exp_date;
        end
        obj.exp_main_dir=pwd;
        obj.exp_dir=fullfile(sup_proj_dir,'Experiments');
        if ~isempty(exp_name)
            obj.exp_save_dir=fullfile(sup_proj_dir,'Experiments',obj.exp_date,obj.exp_name);
            ExpManager.mkdirs(obj.exp_save_dir);
        end
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
    
    function expdir=mkexpdir(obj,exp_name)
        expdir=fullfile(obj.exp_save_dir,exp_name);
        ExpManager.mkdirs(expdir);
    end

    
%     function save_data(obj,data_name,save_name)
%         if nargin<3
%             save_name=data_name;
%         end
%         savePath=fullfile(obj.exp_save_dir,save_name);
%         save savePath data_name;
%         fprintf("save data %s at %s\n",data_name,savePath);    
%     end


%     function load_data(obj,data_name,load_name)
%         if nargin<3
%             load_name=data_name;
%         end
%         load(fullfile(obj.exp_save_dir,load_name),data_name);
%         fprintf("load data %s from %s\n",data_name,loadPath);
%     
%     end
   
end
end