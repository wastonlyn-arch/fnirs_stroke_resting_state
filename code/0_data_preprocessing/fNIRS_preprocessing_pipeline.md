fNIRS 静息态数据预处理自动化流程总结
根据你已完成的操作和经验，下面将整个处理流程梳理成 可复用的标准化步骤。其中部分步骤可以自动化（用 MATLAB 脚本），部分仍需手动操作（如 Homer3 GUI 转换和裁剪）。我会标注出哪些可以整合成脚本，并给出命名规范。

一、项目文件夹结构创建（可脚本化）
目标：按照 BIDS 风格创建规范的文件夹树，包括原始数据、衍生数据、代码等。

已提供的脚本：`create_project_folders.m`（可整合进主脚本）

matlab
```
% create_project_folders.m
% 在项目根目录下创建标准文件夹结构
function create_project_folders(projectRoot)
    if nargin<1, projectRoot = pwd; end
    folders = {'dataset/raw_nirs', 'dataset/raw_snirf', 'derivatives', 'code', 'docs'};
    for i=1:length(folders)
        mkdir(fullfile(projectRoot, folders{i}));
    end
end
```
命名规范：

项目根目录：`fNIRS_stroke_project`

原始数据：`dataset/raw_nirs`（存放原始 .nirs）

转换后数据：`dataset/raw_snirf`（临时存放 .snirf）

被试文件夹：`sub-XX`（如 sub-01, sub-02）

二、文件重命名（可脚本化）
将原始设备导出的 .nirs 文件重命名为 BIDS 兼容格式：
`sub-<ID>_ses-<timepoint>_group-<group>_task-rest.nirs`

命名规范：

被试编号：`sub-01, sub-02 ...`

时间点：`ses-baseline, ses-post`

组别：`group-MT, group-PG, group-MtPg, group-sham`

任务：`task-rest`（静息态）

脚本：`rename_nirs_files.m`（根据原始文件名映射表批量重命名，需要根据你的实际命名规则定制）

三、格式转换：`.nirs → .snirf`（手动 + 半自动）
最佳实践（你已验证）：

将所有 `.nirs` 文件复制到一个临时文件夹（如 D:\temp_nirs）。

打开 MATLAB，运行 `Homer3`。

在 Homer3 中点击 Open，选择该临时文件夹。

Homer3 自动检测所有 .nirs 文件，弹出对话框询问是否转换为 .snirf → 选择 YES（或全选后转换）。

转换完成后，临时文件夹下每个 `.nirs` 旁会生成同名的 `.snirf` 文件。

这一步无法完全自动化（Homer3 的 Nirs2Snirf 函数可能依赖 GUI），但手动操作一次即可。

四、时长检查与长文件准备（可脚本化）
目标：批量读取 .snirf 时长，标记异常（<350 秒 或 >375 秒），并为长文件（>375 秒）创建独立文件夹，供 Homer3 裁剪。

整合脚本：`check_duration_and_prepare_longfiles.m`

matlab
```
% check_duration_and_prepare_longfiles.m
% 功能：
%   1. 读取临时文件夹下所有 .snirf 的时长
%   2. 输出异常文件列表（短于350s或长于375s）
%   3. 为每个时长>375s的文件创建独立子文件夹（便于Segment Tool）
%   4. 将原 .snirf 复制到对应子文件夹中

clear; clc;

% ========== 设置路径 ==========
srcDir = 'D:\fNIRS_mirror_therapy\dataset\.nirs2.snirf';   % 临时文件夹
outDir = fullfile(srcDir, 'to_segment_375plus');           % 长文件准备文件夹

if ~exist(outDir,'dir'), mkdir(outDir); end

snirfFiles = dir(fullfile(srcDir, '*.snirf'));
fprintf('共找到 %d 个 .snirf 文件\n', length(snirfFiles));

results = {};
longFiles = {};

for i = 1:length(snirfFiles)
    filePath = fullfile(snirfFiles(i).folder, snirfFiles(i).name);
    try
        snirf = SnirfClass(filePath);
        snirf.Load();
        t = snirf.data(1).time;
        dur = t(end) - t(1);
        [~,name,~] = fileparts(filePath);
        % 提取被试编号
        subID = regexp(name, 'sub-\d+', 'match');
        if isempty(subID), subID={'unknown'}; else subID=subID{1}; end
        
        results{i,1}=filePath; results{i,2}=dur; results{i,3}=subID;
        
        if dur<350 || dur>375
            fprintf('⚠️ 异常时长 [%.1f s] %s (%s)\n', dur, name, subID{1});
        else
            fprintf('[%3d] %.1f s - %s\n', i, dur, name);
        end
        
        if dur > 375
            longFiles{end+1} = {filePath, dur, subID{1}, name};
            fprintf('    → 将为此长文件创建裁剪文件夹\n');
        end
    catch ME
        warning('读取失败: %s', filePath);
    end
end

% 为长文件建立子文件夹并复制
for j = 1:length(longFiles)
    [srcPath, dur, subID, baseName] = longFiles{j}{:};
    destFolder = fullfile(outDir, baseName);
    if ~exist(destFolder,'dir'), mkdir(destFolder); end
    destPath = fullfile(destFolder, [baseName, '.snirf']);
    if ~exist(destPath,'file')
        copyfile(srcPath, destPath);
        fprintf('已复制: %s -> %s\n', srcPath, destPath);
    end
end

fprintf('\n已完成。请用 Homer3 打开以下文件夹并使用 Segment Tool 裁剪:\n  %s\n', outDir);
```
五、手动裁剪长文件（Homer3 GUI）
在 MATLAB 中运行 Homer3。

点击 Open，选择上一步生成的 `to_segment_375plus` 文件夹。

依次进入每个子文件夹（如 `sub-01_ses-baseline_group-MtPg_task-rest`），加载 `.snirf 文件`。

点击菜单 `Tools` → `Segment SNIRF File`，输入裁剪区间（例如` [0 371]`），确定。

生成的新文件默认名为 原文件名`_seg_1.snirf`。

重复所有长文件的裁剪。裁剪完成后，得到一批` *_seg_1.snirf` 文件，位于各自的子文件夹中。

六、替换裁剪后文件（可脚本化）
将裁剪生成的 *_seg_1.snirf 复制到临时文件夹，覆盖原来的长文件（或替换原临时文件）。

脚本：`replace_segmented_files.m`

matlab
```
% replace_segmented_files.m
% 将 to_segment_375plus 下的所有 *_seg_1.snirf 替换回临时文件夹中的原始文件

clear; clc;
srcDir = 'D:\fNIRS_mirror_therapy\dataset\.nirs2.snirf';
segRoot = fullfile(srcDir, 'to_segment_375plus');

segFiles = dir(fullfile(segRoot, '**', '*_seg_1.snirf'));
for i = 1:length(segFiles)
    segPath = fullfile(segFiles(i).folder, segFiles(i).name);
    [~, segName, ~] = fileparts(segPath);
    origName = strrep(segName, '_seg_1', '');
    origPath = fullfile(srcDir, [origName, '.snirf']);
    if exist(origPath,'file')
        copyfile(segPath, origPath, 'f');
        fprintf('已替换: %s\n', origName);
    else
        warning('原始文件不存在: %s', origPath);
    end
end
```
七、按被试分类放回结构化文件夹（可脚本化）
将临时文件夹中所有最终的 .snirf 文件移动到 dataset/raw_snirf/sub-XX/ 下（或直接放入 raw_fnirs 对应被试文件夹）。

脚本：`organize_snirf_by_subject.m`

matlab
```
% organize_snirf_by_subject.m
clear; clc;
srcDir = 'D:\fNIRS_mirror_therapy\dataset\.nirs2.snirf';      % 最终snirf临时目录
destRoot = 'D:\fNIRS_mirror_therapy\dataset\raw_snirf';       % 目标根目录

snirfFiles = dir(fullfile(srcDir, '*.snirf'));
for i = 1:length(snirfFiles)
    [~,name,~] = fileparts(snirfFiles(i).name);
    subID = regexp(name, 'sub-\d+', 'match');
    if isempty(subID), continue; end
    subDir = fullfile(destRoot, subID{1});
    if ~exist(subDir,'dir'), mkdir(subDir); end
    movefile(fullfile(srcDir, snirfFiles(i).name), fullfile(subDir, [name,'.snirf']));
    fprintf('已移动: %s\n', name);
end
```
八、整合主脚本（建议）
将上述可脚本化的步骤（1、2、4、6、7）整合成一个主脚本 master_pipeline.m，并给出清晰的交互提示，避免重复运行。

示例框架：

matlab
```
% master_pipeline.m
% 完整流程：创建文件夹 → 重命名 → 时长检查 → 准备长文件 → 替换 → 整理

clear; clc;
projectRoot = 'D:\fNIRS_mirror_therapy';

% 1. 创建文件夹结构
create_project_folders(projectRoot);

% 2. 重命名（需要你根据原始文件名自定义映射，此处略）
% rename_nirs_files();

% 3. 格式转换：手动步骤，请用户完成后按任意键继续
disp('请手动用 Homer3 将 .nirs 转换为 .snirf（详见步骤三）');
input('完成后按 Enter 继续...');

% 4. 时长检查并准备长文件
check_duration_and_prepare_longfiles();

% 5. 手动裁剪长文件
disp('请手动用 Homer3 Segment Tool 裁剪长文件（步骤五）');
input('裁剪完成后按 Enter 继续...');

% 6. 替换裁剪后的文件
replace_segmented_files();

% 7. 整理回结构化文件夹
organize_snirf_by_subject();

fprintf('所有处理完成！数据已整理至 %s\n', 'dataset/raw_snirf');
```
总结
| 步骤 | 自动化程度 | 脚本/方法 | 备注 |
| :--- | :--- | :--- | :--- |
| 创建文件夹| 全自动 | `create_project_folders.m` | 可整合进主脚本 |
| 重命名 | 半自动 | `rename_nirs_files.m` | 需根据原始文件名映射表定制 |
| `.nirs` → `.snirf` | 手动（GUI） | Homer3 批量转换 | 一次性操作，无需脚本 |
| 时长检查 & 长文件准备 | 全自动 | `check_duration_and_prepare_longfiles.m` | 自动标记并创建裁剪子文件夹 |
| 裁剪长文件 | 手动（GUI） | Homer3 Segment Tool | 需要人工设置区间（如 [0 371]） |
| 替换裁剪后文件 | 全自动 | `replace_segmented_files.m` | 覆盖原始临时文件 |
| 按被试整理 | 全自动 | `organize_snirf_by_subject.m` | 移动到规范子文件夹 |

将以上脚本保存到项目 `code/` 文件夹中，并根据实际路径修改。以后处理新批次数据时，仅需调整 `projectRoot` 和文件映射表，即可半自动化完成整个预处理流程。
