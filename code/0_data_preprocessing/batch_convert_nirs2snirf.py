#==不适用慧创采集的.nirs文件，MNE 的 read_raw_nirx 函数期望传入的是一个文件夹路径（对应于 NIRx 设备的多文件格式）==
import os
import glob
import mne
#==change root_dir to the directory where your .nirs files are located
root_dir =r'D:\fNIRS_mirror_therapy\dataset\raw_fnirs'
#dir all .nirs files in the root_dir
nirs_files = glob.glob(os.path.join(root_dir,'**','*.nirs'))
print(f'Found {len(nirs_files)} 个.nirs文件')
for nirs_path in nirs_files:
    #生成对应的.snirf文件路径
    snirf_path = os.path.splitext(nirs_path)[0] + '.snirf'
    if os.path.exists(snirf_path):
        print(f'{snirf_path} 已存在，跳过转换')
        continue
    try:
       raw=mne.io.read_raw_nirx(nirs_path, preload=True)
       mne.io.write_raw_snirf(raw, snirf_path)
       print(f'成功转换 {nirs_path} to {snirf_path}')
    except Exception as e:
         print(f'转换 {nirs_path} ')
         print(f'Error: {e}')
print('批量转换完成！')