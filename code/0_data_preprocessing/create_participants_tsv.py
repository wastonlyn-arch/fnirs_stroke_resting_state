import pandas as pd
import os

# 原始分组数据（从你的描述整理）
raw_data = {
    "A": "08,09,14,26,27,30,31,32",   # 号码之间用逗号分隔（原有点和逗号混用，统一为逗号）
    "B": "02,07,11,12,13,18,20",
    "C": "05,06,10,15,16,19,23,24,25",
    "D": "01,03,04,17,21,22,28,29",
}

# 组别映射（请根据实际研究调整）
group_map = {
    "A": "sham",
    "B": "MT",
    "C": "PG",
    "D": "MT&PG",
}

# 输出路径
output_dir = r"D:\fNIRS_mirror therapy\dataset"
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, "participants.tsv")

# 收集所有被试编号和组别
records = []
for group_letter, numbers_str in raw_data.items():
    # 按逗号分割，并去除可能的空格和空值
    numbers = [num.strip() for num in numbers_str.split(",") if num.strip()]
    for num in numbers:
        # 统一格式为两位数字字符串（如 '08'）
        sub_id = f"sub-{int(num):02d}"
        records.append({
            "participant_id": sub_id,
            "name"
            "group": group_map[group_letter],
            "age": "",          # 待补充
            "gender": "",       # 待补充
            "stroke_side": "",  # 待补充
            "pre_FMA": "",      # 待补充
            "post_FMA": "",     # 待补充
            "pre_ARAT": "",     # 待补充
            "post_ARAT": ""     # 待补充
        })

# 创建 DataFrame
df = pd.DataFrame(records)

# 按 participant_id 排序（例如 sub-01, sub-02 ...）
df.sort_values("participant_id", inplace=True)

# 保存为 TSV 文件
df.to_csv(output_file, sep="\t", index=False, encoding="utf-8")

print(f"文件已生成：{output_file}")
print(f"共 {len(df)} 个被试。")
print("\n前5行预览：")
print(df.head().to_string(index=False))