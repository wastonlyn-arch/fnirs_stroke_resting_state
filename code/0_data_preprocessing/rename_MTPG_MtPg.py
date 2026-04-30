import os
import re


def rename_files(directory, dry_run=True):
    """
    递归遍历目录，将文件名中的 'MT&PG' 替换为 'MtPg'
    dry_run=True 时只打印预览，不实际重命名
    """
    for root, dirs, files in os.walk(directory):
        for filename in files:
            if "MT&PG" in filename:
                new_filename = filename.replace("MT&PG", "MtPg")
                old_path = os.path.join(root, filename)
                new_path = os.path.join(root, new_filename)

                if dry_run:
                    print(f"[Dry Run] Would rename:\n  {old_path}\n  -> {new_path}\n")
                else:
                    # 避免覆盖已存在的文件
                    if os.path.exists(new_path):
                        print(f"⚠️ 跳过 {old_path}，因为 {new_path} 已存在")
                    else:
                        os.rename(old_path, new_path)
                        print(f"✓ 已重命名: {filename} -> {new_filename}")


if __name__ == "__main__":
    # 请修改为你的实际根目录（例如包含所有被试文件夹的父目录）
    root_directory = r"D:\fNIRS_mirror_therapy\dataset\raw_fnirs"

    # 先预览（dry_run=True），确认无误后改为 False 执行实际重命名
    rename_files(root_directory, dry_run=False)

    # 确认无误后，将上一行改为：
    # rename_files(root_directory, dry_run=False)