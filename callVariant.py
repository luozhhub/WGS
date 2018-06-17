#!/usr/bin/python
from subprocess import *
import sys
import os
import argparse
from pyflow import WorkflowRunner


"""
useage:
python ./callVariant.py -r csi.chromosome.fa -f1 Flame_1_clean.fq.gz -f2 Flame_2_clean.fq.gz -t /data8/kll/fq/tmp -w /data8/kll/fq -o Flame
"""

class callSnv():
    def __init__(self):
        self.samtools = "/data7/huangyue/software/samtools-1.4.1/samtools"
        self.gatk = "/data8/kll/fq/gatk-package-4.0.4.0-local.jar"
        self.picard = "/data8/kll/fq/picard.jar"
        self.bwa = "/data1/cggi/local/software/bwa-0.7.5a/bwa"
        self.bcftools = "/data7/huangyue/software/bcftools-1.8/bin/bcftools"


    def run(self, cmd, wkdir=None):
        sys.stderr.write("Running %s ...\n" % cmd)
        p = Popen(cmd, shell=True, cwd=wkdir)
        p.wait()
        return p.returncode

    def bwa_index(self, ref=None):
        cmd="%s index %s"%(self.bwa, ref)
        return cmd

    def samtools_index(self, ref=None):
        cmd="%s faidx %s"%(self.samtools, ref)
        return cmd

    def GATK_dict(self, ref=None):
        name_list = ref.split(".")
        name_list = name_list[0:-1]
        output_dict = ".".join(name_list) + ".dict"
        cmd="java -jar %s CreateSequenceDictionary R=%s O=%s"%(self.picard, ref, output_dict)
        return cmd

    def run_bwa(self, core=30, ref=None, fastq_1=None, fastq_2=None, outsam=None, rg=None):
        cmd = "%s mem -R '@RG\tID:1\tPL:ILLUMINA\tSM:%s' -t %s %s %s %s > %s"%(self.bwa, rg, core, ref, fastq_1, fastq_2, outsam)
        return cmd

    def sam_to_bam(self, samfile, bamfile):
        cmd = "%s view -b -S %s -o %s"%(self.samtools, samfile, bamfile)
        return cmd

    def bamsort(self, core=10, mem="3G", sortbamfile=None, bamfile=None, tmpdir=None ,sort_by_name=False):
        if sort_by_name is True:
            cmd="%s sort -n -@ %s -m %s -T %s -o %s %s"%(self.samtools, core, mem, tmpdir, sortbamfile, bamfile)
        else:
            cmd="%s sort -@ %s -m %s -T %s -o %s %s"%(self.samtools, core, mem, tmpdir, sortbamfile, bamfile)
        return cmd

    def bamsort_picard(self, tmpdir=None, sortbamfile=None, bamfile=None):
        cmd = "java -Djava.io.tmpdir=%s -Xms256m -Xmx2048m -jar %s SortSam I=%s O=%s SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT"\
              %(tmpdir, self.picard, bamfile, sortbamfile)
        return cmd

    def sortbamindex(self, sortbam):
        cmd = "%s index %s"%(self.samtools, sortbam)
        return cmd

    #def GATK_addgroup(self):

    def Markduplicate(self, tmpdir=None, input_bam=None, marked_duplicated_bam=None, marked_dup_metrics_txt=None):
        cmd="java -Djava.io.tmpdir=%s -Xms256m -Xmx2048m -jar %s MarkDuplicates I=%s O=%s M=%s VALIDATION_STRINGENCY=LENIENT"\
            %(tmpdir, self.picard, input_bam, marked_duplicated_bam, marked_dup_metrics_txt)
        return cmd

    def Mutect2(self, ref=None, inputbam=None, rg=None, vcffile=None):
        cmd="java -jar %s Mutect2 -R %s -I %s -tumor %s -O %s"\
            %(self.gatk, ref, inputbam, rg, vcffile)
        return cmd

    def select_snv(self, ref=None, input_vcf=None, output_vcf=None):
        cmd="java -jar %s SelectVariants -R %s -V %s --select-type-to-include SNP -O %s"\
            %(self.gatk, ref, input_vcf, output_vcf)
        return cmd

    def mpileup(self, ref=None, input_bam=None, output_bcf=None):
        cmd="%s mpileup -Ob -I -f %s -o %s %s"%(self.bcftools, ref, output_bcf, input_bam)
        return cmd

    def bcftovcf(self, bcfFile=None, vcfFile=None):
        cmd = "%s call -vmO z -o %s %s"%(self.bcftools, vcfFile, bcfFile)
        return cmd


class GATK_workflow(WorkflowRunner, callSnv):
    def __init__(self, fq1=None, fq2=None, tmpdir=None, ref=None, workdir=None, resultdir=None, step=None):
        """
        befor run this work flow, you should make sure you create index for your ref sequence
        bwa index ref
        samtools faidx ref
        java -jar %s CreateSequenceDictionary R=%s O=./GCA_000001405.27_GRCh38.p12_genomic.dict
        :param fq1:
        :param fq2:
        :param tmpdir:
        :param ref:
        """
        callSnv.__init__(self)
        self.workdir = workdir

        self.output_dir = os.path.join(self.workdir, resultdir)
        if os.path.isfile(self.output_dir):
            sys.stderr.write("Output dir %s is a file. Remove it.\n" % self.output_dir)
            os.remove(self.output_dir)
            os.mkdir(self.output_dir)
        elif os.path.isdir(self.output_dir):
            sys.stderr.write("Output dir %s exists.\n" % self.output_dir)
        else:
            os.mkdir(self.output_dir)
        self.fq1 = os.path.join(self.workdir, fq1)
        self.fq2 = os.path.join(self.workdir, fq2)
        self.sample_id = fq1.split("_")[0]
        self.sam = os.path.join(self.output_dir, self.sample_id + ".sam")
        self.bam = os.path.join(self.output_dir, self.sample_id + ".bam")
        self.sortbam = os.path.join(self.output_dir, self.sample_id + ".sort.bam")
        self.mkdupsortbam = os.path.join(self.output_dir, self.sample_id + ".mkdup.sort.bam")
        self.mkmatrix = os.path.join(self.output_dir, self.sample_id + ".matrix.txt")
        self.bai = self.mkdupsortbam + ".bai"
        self.tmpbasedir = tmpdir
        self.tmpdir = os.path.join(tmpdir, self.sample_id)
        self.bcffile = os.path.join(self.output_dir, self.sample_id + ".bcf")
        self.vcffile = os.path.join(self.output_dir, self.sample_id + ".vcf.gz")
        self.snv_vcffile = os.path.join(self.output_dir, self.sample_id + ".snv.vcf.gz")
        self.ref = os.path.join(self.workdir, ref)
        self.step = step


    def workflow(self):
        #step one
        if self.step < 1:
            MapFastqToSam = self.addTask("MapFastqToSam", self.run_bwa(core=10, ref=self.ref, fastq_1=self.fq1, fastq_2=self.fq2, outsam=self.sam, rg=self.sample_id))
        else:
            MapFastqToSam = self.addTask("MapFastqToSam")

        #step_two
        if self.step < 2:
            SamToBam = self.addTask("SamToBam", self.sam_to_bam(samfile=self.sam, bamfile=self.bam), dependencies=[MapFastqToSam])
        else:
            SamToBam = self.addTask("SamToBam", dependencies=[MapFastqToSam])

        #step three
        if self.step < 3:
            SortBam = self.addTask("SortBam", self.bamsort_picard(tmpdir=self.tmpbasedir, sortbamfile=self.sortbam, bamfile=self.bam), dependencies=[SamToBam])
        else:
            SortBam = self.addTask("SortBam", dependencies=[SamToBam])

        #step four
        if self.step < 4:
            markDuplicate = self.addTask("markduplicate", self.Markduplicate(tmpdir=self.tmpbasedir, input_bam=self.sortbam, marked_duplicated_bam=self.mkdupsortbam, marked_dup_metrics_txt=self.mkmatrix), dependencies=[SortBam])
        else:
            markDuplicate = self.addTask("markduplicate", dependencies=[SortBam])

        # step_five
        if self.step < 5:
            BamIndex = self.addTask("BamIndex", self.sortbamindex(self.mkdupsortbam), dependencies=[markDuplicate])
        else:
            BamIndex = self.addTask("BamIndex", dependencies=[markDuplicate])

        #step_six
        if self.step < 6:
            mpileup = self.addTask("mpileup", self.mpileup(ref=self.ref, input_bam=self.mkdupsortbam, output_bcf=self.bcffile), dependencies=[BamIndex])
        else:
            mpileup = self.addTask("mpileup", dependencies=[BamIndex])

        #step seven
        bcftovcf = self.addTask("bcftovcf", self.bcftovcf(bcfFile=self.bcffile, vcfFile=self.vcffile), dependencies=[mpileup])
        """
        if self.step < 6:
            mutect2 = self.addTask("mutect2", self.Mutect2(ref=self.ref, inputbam=self.mkdupsortbam, rg=self.sample_id, vcffile=self.vcffile), dependencies=[BamIndex])
        else:
            mutect2 = self.addTask("mutect2", dependencies=[BamIndex])
        """



        #step seven
        selectVariant = self.addTask("selectVariant", self.select_snv(ref=self.ref, input_vcf=self.vcffile, output_vcf=self.snv_vcffile), dependencies=[mutect2])
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument("-c", "--core", type=int, default=10)
    ap.add_argument("-r", "--ref", type=str, required=True)
    ap.add_argument("-f1", "--fastq_1", type=str, required=True)
    ap.add_argument("-f2", "--fastq_2", type=str, required=True)
    ap.add_argument("-t", "--tmp", type=str, required=True)
    ap.add_argument("-w", "--wkdir", type=str, required=True)
    ap.add_argument("-o", "--output_dir", type=str, required=True)
    ap.add_argument("-s", "--step", type=int, default=0,
                    help='0: start from very begining\n' \
                         '1: convert the sam to bam file\n' \
                         '2: sort bam file\n' \
                         '3: markduplicates\n' \
                         '4: make bam index file .bai\n'\
                         '5: use mutect2 call snv'
                    )
    args = ap.parse_args()
    wf = GATK_workflow(fq1=args.fastq_1, fq2=args.fastq_2, tmpdir=args.tmp, ref=args.ref, workdir=args.wkdir, resultdir=args.output_dir, step=args.step)
    retval = wf.run(mode="local", nCores=10)
    sys.exit(retval)


