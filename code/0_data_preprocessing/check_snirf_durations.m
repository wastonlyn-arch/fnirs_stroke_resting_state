% check_and_prepare_for_segment.m
% 功能：
%   1. 批量读取 .snirf 文件，获取每个文件的时长（秒）
%   2. 筛选出时长 <350 或 >375 的文件，并打印异常汇总
%   3. 对于时长 >375 的文件，自动创建文件夹并复制文件，以便在 Homer3 中统一裁剪

clear; clc;

% ========== 请修改为你的 .snirf 文件所在根目录 ==========
rootDir = 'D:\fNIRS_mirror_therapy\dataset\.nirs2.snirf';

% 定义输出目录，用于存放需要裁剪的文件（将按原文件名建立子文件夹）
outputDir = fullfile(rootDir, 'to_segment_375plus');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% 递归查找所有 .snirf 文件
snirfFiles = dir(fullfile(rootDir, '**', '*.snirf'));

if isempty(snirfFiles)
    error('未找到任何 .snirf 文件');
end

fprintf('共找到 %d 个 .snirf 文件\n\n', length(snirfFiles));

% 存储结果和异常文件信息
results = {};
abnormalFiles = {};   % 存储异常文件的信息 {path, duration, subID}
longFiles = {};       % 存储时长 >375 的文件信息，用于后续准备文件夹

for i = 1:length(snirfFiles)
    filePath = fullfile(snirfFiles(i).folder, snirfFiles(i).name);
    try
        % 读取 .snirf 文件
        snirf = SnirfClass(filePath);
        t = snirf.data(1).time;
        duration = t(end) - t(1);
        
        % 提取被试编号（假设文件名中包含 sub-XX 格式）
        [~, name, ~] = fileparts(filePath);
        subID = regexp(name, 'sub-\d+', 'match');
        if isempty(subID)
            subID = {'unknown'};
        else
            subID = subID{1};
        end
        
        results{i, 1} = filePath;
        results{i, 2} = duration;
        results{i, 3} = subID;
        
        % 筛选异常（<350 或 >375）
        if duration < 350 || duration > 375
            abnormalFiles{end+1} = {filePath, duration, subID};
            fprintf('⚠️ 异常时长 [%6.1f 秒]  %s  (被试: %s)\n', duration, filePath, subID);
        else
            fprintf('[%3d] %6.1f 秒  -  %s\n', i, duration, filePath);
        end
        
        % 单独记录长文件（>375）以便后续准备文件夹
        if duration > 375
            longFiles{end+1} = {filePath, duration, subID, name};  % name 不含扩展名
            fprintf('     → 该文件时长 >375 秒，将为其创建裁剪文件夹\n');
        end
        
    catch ME
        warning('读取失败: %s\n错误: %s', filePath, ME.message);
        results{i, 1} = filePath;
        results{i, 2} = NaN;
        results{i, 3} = 'unknown';
    end
end

% ========== 为时长 >375 的文件创建文件夹并复制文件 ==========
if ~isempty(longFiles)
    fprintf('\n========== 为时长 >375 的文件准备裁剪文件夹 ==========\n');
    for j = 1:length(longFiles)
        srcPath = longFiles{j}{1};
        fileName = longFiles{j}{4};
        destFolder = fullfile(outputDir, fileName);
        if ~exist(destFolder, 'dir')
            mkdir(destFolder);
        end
        % 复制文件（将 .snirf 复制到子文件夹中，保持原文件名）
        [~, fname, ext] = fileparts(srcPath);
        destPath = fullfile(destFolder, [fname, ext]);
        if ~exist(destPath, 'file')
            copyfile(srcPath, destPath);
            fprintf('已复制: %s -> %s\n', srcPath, destPath);
        else
            fprintf('文件已存在，跳过复制: %s\n', destPath);
        end
    end
    fprintf('所有长文件已准备就绪，请用 Homer3 打开以下文件夹进行裁剪：\n  %s\n', outputDir);
    fprintf('裁剪目标：统一长度至 371 秒（可在 Segment Tool 中设置 Start=0, End=371）\n');
else
    fprintf('\n没有发现时长 >375 秒的文件，无需准备裁剪文件夹。\n');
end

% ========== 输出异常汇总 ==========
fprintf('\n========== 异常时长汇总（<350 或 >375） ==========\n');
if isempty(abnormalFiles)
    fprintf('未发现异常时长的文件。\n');
else
    fprintf('共发现 %d 个异常文件：\n', length(abnormalFiles));
    for i = 1:length(abnormalFiles)
        fprintf('  %s  (%.1f 秒)  被试：%s\n', abnormalFiles{i}{1}, abnormalFiles{i}{2}, abnormalFiles{i}{3});
    end
end

% ========== 保存详细结果到 CSV ==========
outputTable = table(results(:,1), cell2mat(results(:,2)), results(:,3), ...
    'VariableNames', {'FilePath', 'DurationSec', 'SubjectID'});
writetable(outputTable, 'snirf_durations_all.csv');
fprintf('\n所有文件时长已保存至 snirf_durations_all.csv\n');