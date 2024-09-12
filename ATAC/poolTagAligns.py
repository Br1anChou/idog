#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @File       :poolTagAligns.py
# @Time       :2024/9/1 14:54
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :atac_pipeline
# @Version    :python 3.10.6
# @Usage      :python poolTagAligns.py <output_name> <input_file1> <input_file2> ...

import os
import subprocess
import sys
from sys import exit

def pool_ta(output_name, tas, out_dir='.'):
	if len(tas) > 1:
		# 设置输出文件的完整路径
		pooled_ta: str = os.path.join(out_dir, f'{output_name}.pooled.tn5.tagAlign.gz')

		# 构建命令
		cmd = f'zcat -f {" ".join(tas)} ｜ gzip -nc > {pooled_ta}'

		# 执行命令
		subprocess.run(cmd, shell=True, check=True)

		print(f'Pooled file created: {pooled_ta}')
		return pooled_ta
	else:
		print('Needs at least two TAs (or BEDs) to be pooled.')
		return None

if __name__ == '__main__':
	# 确保至少有一个输出文件名和两个输入文件名
	if len(sys.argv) < 3:
		print("Usage: python3 poolTagAligns <output_name> <input_file1> <input_file2> ...")
		exit(1)
	else:
		pass

	# 第一个参数是输出文件名
	output_name = sys.argv[1]

	# 后续参数是输入文件列表
	tas = sys.argv[2:]

	# 调用合并函数
	pool_ta(output_name=output_name, tas=tas, out_dir='./')
