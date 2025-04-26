clc;clear all
kpath=load('BTE.qpoints');
infile_qpoints=[kpath(1:end-1,4:end),kpath(2:end,4:end)];

% 指定文件名
filename = 'infile.qpoints_dispersion';

% 打开文件进行写入操作
fileID = fopen(filename, 'w');

% 检查文件是否成功打开
if fileID == -1
    error('无法打开文件进行写入。');
end

% 将文本写入文件
fprintf(fileID, '%s\n', 'CUSTOM');
% 写入第二行
fprintf(fileID, '%d\n', 2);

fprintf(fileID, '%d\n', size(infile_qpoints,1));

for i=1:size(infile_qpoints,1)
    fprintf(fileID, '%d %d %d %d %d %d ', infile_qpoints(i,:));
    fprintf(fileID, '%s\n', ' X X');
end

% 关闭文件
fclose(fileID);  