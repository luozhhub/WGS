#!/usr/bin/python
from subprocess import *
import sys
import os
import argparse
import gzip
import re


class TRANS():
    def __init__(self, outfile=None):
        self.inputFile = "combime_file.snp"
        self.sample_list = ["19P", "Flame", "WBY", "m2", "ZK", "YLK", "SO3", "JY", "HKC"]
        self.header = "\t9\t"
        self.outputHandle = open(outfile, "w")



    def write_line(self, inputFile=None, sample=None, leng=None, outHandle=None, filter=False):
        prefix = sample
        while len(prefix) < leng:
            prefix = prefix + " "
        line_text = prefix
        fileHandle = open(inputFile, "r")
        index_n = self.sample_list.index(sample)
        nlen = 0
        for line in fileHandle:
            chars = line.strip("\n").split("\t")[2]
            if filter is True:
                array = list(chars)
                if len(set(array)) == 1:
                    continue
            nlen = nlen + 1
            char = chars[index_n]
            line_text = line_text + char
        if self.header != "":
            self.header = self.header + str(nlen) + "\n"
            outHandle.write(self.header)
            self.header = ""
        outHandle.write(line_text + "\n")
        fileHandle.close()

    def run(self, filter=False):
        for sam in self.sample_list:
            self.write_line(inputFile=self.inputFile, sample=sam, leng=22, outHandle=self.outputHandle, filter=filter)
        self.outputHandle.close()


if __name__ == '__main__':
    trans = TRANS(outfile="transform_all.snp")
    trans.run(filter=False)

    trans_1 = TRANS(outfile="transform_notall.snp")
    trans_1.run(filter=True)

