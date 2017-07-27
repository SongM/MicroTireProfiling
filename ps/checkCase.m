function mat_file = checkCase(data_folder, case_name)
    %% check if case is active
    mat_file = [data_folder, '\', case_name, '\', case_name, '.mat'];
    if (exist(mat_file,'file'))
        % if exist, display current situation
        load(mat_file);
        fprintf(1,['    function_folder = \n']);    disp(function_folder);
        fprintf(1,['    case_name = \n']);  disp(case_name);
        job_done
        convert_task
    else
        % if not exist, create
        fprintf(['case ', case_name, ' does not exist, create the case for it.\n']);
        function_folder = 'C:\Users\s\Google Drive\projects\honda_20170118\functions';
        folder_path = [data_folder, '\', case_name];
        mkdir([folder_path,'\processed_data']);
%         mkdir([folder_path,'\processed_data\back_up']);
        figure_count = 100;
        convert_task = [];
        job_done = {'0.just start'};
        param = [];
        error = [];
        individual_mat_file = [];
        
        saveMatFile();
    end
end

