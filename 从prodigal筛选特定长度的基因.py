import os
import re
import argparse

parser = argparse.ArgumentParser(description="You should add those parameters")
parser.add_argument("-i","--input", type = str, help = "The input fasta file")
parser.add_argument("-min","--min_length", type = int, help = "The minimal length of sequences")
parser.add_argument("-max","--max_length", type = int, default = 0, help = "The maximal length of sequences")
parser.add_argument("-o","--output", type = str, help = "The output fasta file")
#args = parser.parse_args()
args =parser.parse_known_args()[0]

# 输入的fasta文件转换成字典
file_open = open(args.input,"r")
file_read = file_open.readlines()
res_dict = {}
for line in file_read:
    if re.match(">",line):
        #print(line)
        res_dict[line] = ""
        flag = line
    else:
        res_dict[flag] = res_dict[flag] + line

# 根据长度进行筛选
# 根据长度进行筛选
final_results = []
for key, value in res_dict.items():
    if args.max_length == 0:
        if len(value) >= args.min_length:
            final_results.append(key)
            final_results.append(value)
    else:
        if args.min_length <= len(value) <= args.max_length:
            final_results.append(key)
            final_results.append(value)

# 保存文件
file_new = open(args.output,"w")
for i in final_results:
    file_new.write(i + "\n")

file_new.close()