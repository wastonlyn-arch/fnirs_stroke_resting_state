import os
import shutil

# 原始数据根目录（包含 con1/, con2/ 等文件夹）
source_root = r"D:\fNIRS_mirror therapy"

# 目标根目录（你的项目 dataset/raw_fNIRS 文件夹）
target_root = r"D:\fNIRS_mirror therapy\dataset\raw_fnirs"

# 定义组别映射（请根据你的实际情况修改）
group_map = {
    "A": "sham",    # 假设 A -> 镜像疗法组
    "B": "MT",    # 假设 B -> 手套组
    "C": "PG",  # 假设 C -> 联合干预组
    "D": "MT&PG",  # 假设 D -> 对照组
}

# 时间点映射（当前只有 con1 和 con2）
timepoint_map = {
    "con1": "baseline",
    "con2": "post",   # 如果你有干预后数据
}

# 遍历所有 .nirs 文件
for root, dirs, files in os.walk(source_root):
    for file in files:
        if not file.endswith(".nirs"):
            continue
        # 解析路径: 期望格式为 .../conX/组别/数字.nirs
        parts = root.split(os.sep)
        # 假设倒数第二层是组别（如 "A"），倒数第三层是时间点（如 "con1"）
        if len(parts) < 2:
            print(f"跳过无法解析路径: {root}")
            continue
        timepoint = parts[-2]   # 例如 "con1"
        group_letter = parts[-1]   # 例如 "A"
        subj_num = file.replace(".nirs", "")   # 例如 "08"

        # 跳过非数字的被试编号
        if not subj_num.isdigit():
            print(f"跳过非数字编号: {file}")
            continue

        # 映射组别和时间点
        group_name = group_map.get(group_letter, "unknown")
        tp_name = timepoint_map.get(timepoint, timepoint)

        # 构造新文件名和路径
        new_filename = f"sub-{int(subj_num):02d}_ses-{tp_name}_group-{group_name}_task-rest.nirs"
        target_subdir = os.path.join(target_root, f"sub-{int(subj_num):02d}")
        os.makedirs(target_subdir, exist_ok=True)
        target_path = os.path.join(target_subdir, new_filename)

        # 复制（或移动）文件
        shutil.copy2(os.path.join(root, file), target_path)  # 复制；若想移动则用 shutil.move
        print(f"已处理: {os.path.join(root, file)} -> {target_path}")

print("重命名完成！")