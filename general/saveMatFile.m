save(mat_file);
c = clock;
disp(datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6))));
timestamp = ['_',num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5))];
mat_file_w_timestamp = [folder_path,'\processed_data\',case_name,timestamp,'.mat'];
clear c timestamp
save(mat_file_w_timestamp);
clear mat_file_w_timestamp
clear c;
diary off;
diary on;
