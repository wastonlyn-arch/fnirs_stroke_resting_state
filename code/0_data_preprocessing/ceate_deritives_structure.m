% create_derivatives_structure.m
% 根据 BIDS 规范创建 fNIRS 分析的 derivatives 文件夹结构

clear; clc;

% ========== 用户配置区域 ==========
% 项目根目录（包含 derivatives 文件夹的父目录）
rootPath = 'D:/fNIRS_mirror_therapy';   % 请修改为你的实际路径


% 自动生成连续编号的被试（若使用方式二，请注释上面的 subjectList）
 numSubjects = 32;
 subjectList = arrayfun(@(x) sprintf('sub-%02d', x), 1:numSubjects, 'UniformOutput', false);

% ========== 无需修改以下代码 ==========
derivRoot = fullfile(rootPath, 'derivatives', 'homer');
if ~exist(derivRoot, 'dir')
    mkdir(derivRoot);
    fprintf('创建目录: %s\n', derivRoot);
end

% 创建每个被试的子目录
for i = 1:length(subjectList)
    subDir = fullfile(derivRoot, subjectList{i});
    if ~exist(subDir, 'dir')
        mkdir(subDir);
        fprintf('创建目录: %s\n', subDir);
    end
    
    % 创建功能连接、网络、ROI、质量控制等子目录
    subDirs = {'func_conn', 'network', 'roi', 'qc'};
    for j = 1:length(subDirs)
        subSubDir = fullfile(subDir, subDirs{j});
        if ~exist(subSubDir, 'dir')
            mkdir(subSubDir);
            fprintf('创建目录: %s\n', subSubDir);
        end
    end
end

% 创建空的 .processing_stream_config.cfg 文件（位于 derivatives/homer/ 下）
cfgFile = fullfile(derivRoot, '.processing_stream_config.cfg');
if ~exist(cfgFile, 'file')
    fid = fopen(cfgFile, 'w');
    fclose(fid);
    fprintf('创建配置文件: %s\n', cfgFile);
else
    fprintf('配置文件已存在: %s\n', cfgFile);
end

fprintf('\n衍生文件夹结构创建完成！\n');
fprintf('根目录: %s\n', derivRoot);