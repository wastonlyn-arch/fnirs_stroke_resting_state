import os
base_dir="fNIRS_mirror therapy"
subdirs=[
    "dataset/raw_fnirs",
    "dataset/raw_clinical",
    "derivatives/subject_level_stats",
    "derivatives/group_level_stats",
    "derivatives/results/behavioral_analysis",
    "derivatives/results/fnirs_analysis",
    "derivatives/results/correlation",
    "code/0_data_preprocessing",
    "code/1_fnirs_preprocess",
    "code/2_behavioral_stats",
    "code/3_fnirs_analysis",
    "code/4_correlation_analysis",
    "code/5_figures_tables",
    "code/utils",
    "docs"
]
#循环创建每个子文件夹
for folder in subdirs:
    full_path=os.path.join(base_dir,folder)
    os.makedirs(full_path,exist_ok=True)
    print(f"created:{folder}")
#创建两个占位文件
readme_path=os.path.join(base_dir,"README.md")
with open(readme_path,"w") as f:
    f.write("#卒中手功能康复析因研究\n\n项目说明...")
print(f"created:{readme_path}")

req_path=os.path.join(base_dir,"requirements.txt")
with open(req_path,"w") as f:
    f.write("#依赖包列表\npandas\nnumpy\nmne\nmatplotlib\nseaborn")
print(f"created:{req_path}")