%%
function_folder = 'C:\Users\s\Google Drive\projects\honda_20170118\functions';
addpath(genpath(function_folder));

rmpath(genpath([function_folder,'\','New Folder']));getStarted;

case_name = '20170724';
data_folder = 'C:\Users\s\Google Drive\projects\honda_20170118\laser_photo';

mat_file = checkCase(data_folder, case_name);
load(mat_file);

%%
param.stitching.ref_groove_region_Z =  [-15,5];
param.stitching.tire_surface_region_Z = [-100,50];
param.stitching.groove_region_half_size_Z = 3;
param.stitching.N_groove = 4;
param.stitching.fitting_err_max_th = 0.5;
% param.stitching.
saveMatFile();
seperateGrooveAndSurfacePixel(mat_file);



c = clock;
disp(datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6))));
timestamp = ['_',num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5))];
individual_mat_file.stitching = [data_folder, '\', case_name, '\processed_data\', 'stitching',timestamp,'.mat'];
clear c timestamp;
saveMatFile();



findTireAxis(mat_file);


for i=1:20
    load(['ps_res_',num2str(i),'.mat']);
    ps_set_res(i) = ps_res;
end

stitchPsRes( mat_file, ps_set_res )



load(mat_file);
load(individual_mat_file.stitching,'tire');
final_res = getFinalResults(tire);

c = clock;
disp(datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6))));
timestamp = ['_',num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5))];
individual_mat_file.final_res = [data_folder, '\', case_name, '\processed_data\', 'final_res',timestamp,'.mat'];
clear c timestamp;
save(individual_mat_file.final_res,'final_res');
