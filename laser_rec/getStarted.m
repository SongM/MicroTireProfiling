clc;
warning off;
c = clock;
diary on;
diary(['C:\Users\s\Dropbox\code\matlab\log\log','_',num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'.out']);
disp(datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6))));
disp('Have a good day!');
clear all;close all;
diary;