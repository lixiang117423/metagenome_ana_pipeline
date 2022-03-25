import os
import re
import argparse

parser = argparse.ArgumentParser(description="You should add those parameters")
parser.add_argument("-i","--input", type = str, help = "The input fasta file")
parser.add_argument("-label","--label", type = str, help = "The label used to add on the old ID")
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
file_open.close()

# 添加标签
final_results = []
for key, value in res_dict.items():
    key = key.replace(">", ">" + args.label + "|")
    final_results.append(key)
    final_results.append(value)

# 保存文件
file_new = open(args.output,"w")
for i in final_results:
    file_new.write(i + "\n")

file_new.close()