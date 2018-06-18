#!/usr/bin/python
from subprocess import *
import sys
import os
import argparse
import gzip
import re


def read_gz_file(path):
    if os.path.exists(path):
        with open(path, "r") as pf:
            for line in pf:
                yield line
    else:
        print("the path %s is not exits"%(path))


sample_list = ["19P", "Flame", "WBY", "m2", "ZK", "YLK", "SO3", "JY", "HKC"]
chrs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrUn']
file_list=[]
file_handle = {}
for sample in sample_list:
    file_path_1 = os.path.join("/data8/kll/fq", sample, sample + ".convert.snp")
    file_list.append(file_path_1)
    file_handle[file_path_1] = read_gz_file(file_path_1)

def get_items(file_path=None, chr=None, pos=None):
    con = file_handle[file_path]
    items = con.next().strip("\n").split("\t")
    items[1] = int(items[1])
    while chrs.index(items[0]) < chr:
        items = con.next().strip("\n").split("\t")
        items[1] = int(items[1])
    while items[1] < pos:
        items = con.next().strip("\n").split("\t")
        items[1] = int(items[1])
    #exit(1)
    return items

output = open("combime_file.snp", "w")
store_file_path = ""
store_items = []
logFile = open("combine.snp.log", "w")
logFile_1 = open("combine_1.snp.log", "w")

for line in open(file_list[0], "r"):
    items = line.strip("\n").split("\t")
    chr = chrs.index(items[0])
    pos = int(items[1])
    if store_items != []:
        if pos < store_items[1]:
            continue
    letter = items[2]
    stri = letter
    #logFile_1.write(str(pos) + "\t")
    for file_path in file_list[1:]:
        #logFile_1.write(file_path + "\n")
        if store_file_path == file_path:
            items_x = store_items
            if items_x[1] < pos:
                items_x = get_items(file_path=file_path, chr=chr, pos=pos)
        else:
            items_x = get_items(file_path=file_path, chr=chr, pos=pos)
            logFile_1.write(str(items_x))
            logFile_1.write("\t" +file_path +"\n")
            #print file_path
            #print items_x[1]
        if items_x[1] == pos:
            stri = stri + items_x[2]
        if items_x[1] > pos:
            store_file_path = file_path
            store_items = items_x
            break
    if len(stri) < 9:
        logFile.write(str(store_items) + "\n")
        logFile.write(store_file_path + "\t" + str(store_items[1]) + "\n")
        #exit(1)
        continue
    else:
        output.write(chrs[chr] + "\t" + str(pos) + "\t" + stri + "\n")
        store_file_path = ""
        store_items = []

