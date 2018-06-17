#!/usr/bin/python
from subprocess import *
import sys
import os
import argparse
from pyflow import WorkflowRunner


"""
useage:
python ./bammappy.py -r hs37d5.fa -f1 WGC090555D_combined_R1.fastq.gz -f2 WGC090555D_combined_R2.fastq.gz -t /simm/scratch/tmp -w /simm/Sunset/workflow/scratch/testdir/WGC090555D -o WGC090555D_pdx
"""

class callSnv():
    def __init__(self):
        self.samtools = "/simm/program/bin/samtools"
        self.bamcmp = "/simm/home/luozhihui/bamcmp_map/bamcmp/build/bamcmp"
        self.gatk = "/simm/Sunset/workflow/scratch/testdir/gatk-package-4.0.4.0-local.jar"
        self.picard = "/simm/home/luozhihui/my_source_code/picard.jar"
        self.bwa = "/simm/home/luozhihui/program/bwa-0.7.12/bwa"


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
        cmd="java -jar %s CreateSequenceDictionary R=%s O=./GCA_000001405.27_GRCh38.p12_genomic.dict"%(self.picard, ref, output_dict)
        return cmd

    def run_bwa(self, core=30, ref=None, fastq_1=None, fastq_2=None, outsam=None, rg=None):
        cmd = "%s mem -R '@RG\tID:1\tPL:ILLUMINA\tSM:%s' -t %s %s %s %s > %s"%(self.bwa, rg, core, ref, fastq_1, fastq_2, outsam)
        return cmd

    def sam_to_bam(self, samfile, bamfile):
        cmd = "%s view -b -S %s -o %s"%(self.samtools, samfile, bamfile)
        return cmd

    def bamsort(self, core=10, mem="3G", sortbamfile=None, bamfile=None, tmpdir=None ,sort_by_name=False):
        if sort_by_name is True:
            cmd="%s sort -n -@ %s -m %s -T %s -o %s %s"%(self.samtools, core, mem, sortbamfile, bamfile, tmpdir)
        else:
            cmd="%s sort -@ %s -m %s -T %s -o %s %s"%(self.samtools, core, mem, sortbamfile, bamfile, tmpdir)
        return cmd

    def sortbamindex(self, sortbam):
        cmd = "%s index %s"%(self.samtools, sortbam)
        return cmd

    #def GATK_addgroup(self):

    def Markduplicate(self, tmpdir=None, input_bam=None, marked_duplicated_bam=None, marked_dup_metrics_txt=None):
        cmd="java -Djava.io.tmpdir=%s -jar -Xms1024m -Xmx10240m %s MarkDuplicates I=%s O=%s M=%s VALIDATION_STRINGENCY=LENIENT"\
            %(tmpdir, self.picard, input_bam, marked_duplicated_bam, marked_dup_metrics_txt)
        return cmd

    def Mutect2(self, tmpdir=None, ref=None, inputbam=None, rg=None, vcffile=None):
        cmd="java -Djava.io.tmpdir=%s -jar -Xms1024m -Xmx10240m %s Mutect2 -R %s -I %s -tumor %s -O %s"\
            %(tmpdir, self.gatk, ref, inputbam, rg, vcffile)
        return cmd

    def select_snv(self, tmpdir=None, ref=None, input_vcf=None, output_vcf=None):
        cmd="java -Djava.io.tmpdir=%s -jar -Xms1024m -Xmx10240m %s SelectVariants -R %s -V %s --select-type-to-include SNP -O %s"\
            %(tmpdir, self.gatk, ref, input_vcf, output_vcf)
        return cmd

    def bam_cmp(self, human_sortbamfile=None, mouse_sortbamfile=None, humanonlybam=None, humanbetterbam=None, mouseonlybam=None, mousebetterbam=None, humanlossbam=None, mouselossbam=None):
        cmd = "%s -n -1 %s -2 %s -a %s -A %s -b %s -B %s -C %s -D %s -s mapq"\
              %(self.bamcmp, human_sortbamfile, mouse_sortbamfile, \
                  humanonlybam, humanbetterbam, mouseonlybam, mousebetterbam, humanlossbam, mouselossbam)
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
        self.vcffile = os.path.join(self.output_dir, self.sample_id + ".vcf.gz")
        self.snv_vcffile = os.path.join(self.output_dir, self.sample_id + ".snv.vcf.gz")
        self.ref = os.path.join(self.workdir, ref)
        self.step = step


    def workflow(self):
        #step one
        if self.step < 1:
            MapFastqToSam = self.addTask("MapFastqToSam", self.run_bwa(core=30, ref=self.ref, fastq_1=self.fq1, fastq_2=self.fq2, outsam=self.sam, rg=self.sample_id))
        else:
            MapFastqToSam = self.addTask("MapFastqToSam")

        #step_two
        if self.step < 2:
            SamToBam = self.addTask("SamToBam", self.sam_to_bam(samfile=self.sam, bamfile=self.bam), dependencies=[MapFastqToSam])
        else:
            SamToBam = self.addTask("SamToBam", dependencies=[MapFastqToSam])

        #step three
        if self.step < 3:
            SortBam = self.addTask("SortBam", self.bamsort(core=10, mem="3G", sortbamfile=self.sortbam, bamfile=self.bam, tmpdir=self.tmpdir), dependencies=[SamToBam])
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
            mutect2 = self.addTask("mutect2", self.Mutect2(tmpdir=self.tmpbasedir, ref=self.ref, inputbam=self.mkdupsortbam, rg=self.sample_id, vcffile=self.vcffile), dependencies=[BamIndex])
        else:
            mutect2 = self.addTask("mutect2", dependencies=[BamIndex])

        #step seven
        selectVariant = self.addTask("selectVariant", self.select_snv(ref=self.ref, input_vcf=self.vcffile, output_vcf=self.snv_vcffile), dependencies=[mutect2])

class bamcmp_workflow(GATK_workflow):
    def __init__(self, fq1=None, fq2=None, tmpdir=None, ref=None, workdir=None, resultdir=None, step=None):
        GATK_workflow.__init__(self, fq1=fq1, fq2=fq2, tmpdir=tmpdir, ref=ref, workdir=workdir, resultdir=resultdir, step=step)
        self.mouseSam = os.path.join(self.output_dir, self.sample_id + ".mouse.sam")
        self.mouseBam = os.path.join(self.output_dir, self.sample_id + ".mouse.bam")
        self.mouseref = os.path.join(self.workdir,"GCA_000001635.8_GRCm38.p6_genomic.fna.gz")
        self.humanBam = os.path.join(self.output_dir, self.sample_id + ".human.bam")

        self.sortMouseBam = os.path.join(self.output_dir, self.sample_id + ".sortN.mouse.bam")
        self.sortHumanBam = os.path.join(self.output_dir, self.sample_id + ".sortN.human.bam")
        self.humanOnlyBam = os.path.join(self.output_dir, self.sample_id + ".human.only.bam")
        self.humanBetterBam = os.path.join(self.output_dir, self.sample_id + ".human.better.bam")
        self.mouseOnlyBam = os.path.join(self.output_dir, self.sample_id + ".mouse.only.bam")
        self.mouseBetterBam = os.path.join(self.output_dir, self.sample_id + ".mouse.better.bam")
        self.humanLossBam = os.path.join(self.output_dir, self.sample_id + ".human.loss.bam")
        self.mouseLossBam = os.path.join(self.output_dir, self.sample_id + ".mouse.loss.bam")

        self.humanBetterSortBam = os.path.join(self.output_dir, self.sample_id + ".human.better.sort.bam")


    def workflow(self):
        #step one
        if self.step <1:
            MapFastqToSam = self.addTask("MapFastqToSam",
                                     self.run_bwa(core=30, ref=self.mouseref, fastq_1=self.fq1, fastq_2=self.fq2,
                                                  outsam=self.mouseSam, rg=self.sample_id))
        else:
            MapFastqToSam = self.addTask("MapFastqToSam")

        #step two
        if self.step < 2:
            SamToBam = self.addTask("SamToBam", self.sam_to_bam(samfile=self.mouseSam, bamfile=self.mouseBam),
                                dependencies=[MapFastqToSam])
        else:
            SamToBam = self.addTask("SamToBam", dependencies=[MapFastqToSam])

        #step three
        if self.step < 3:
            SortBamN_1 = self.addTask("SortBamN_1", self.bamsort(core=10, mem="3G", sortbamfile=self.sortMouseBam, bamfile=self.mouseBam,
                                                       tmpdir=self.tmpdir, sort_by_name=True), dependencies=[SamToBam])
            SortBamN_2 = self.addTask("SortBamN_2", self.bamsort(core=10, mem="3G", sortbamfile=self.sortHumanBam, bamfile=self.humanBam,
                                                                 tmpdir=self.tmpdir, sort_by_name=True), dependencies=[SamToBam])
        else:
            SortBamN_1 = self.addTask("SortBamN_1", dependencies=[SamToBam])
            SortBamN_2 = self.addTask("SortBamN_2", dependencies=[SamToBam])



        #step four
        if self.step < 4:
            BamCmp = self.addTask("BamCmp", self.bam_cmp(human_sortbamfile=self.sortHumanBam, mouse_sortbamfile=self.sortMouseBam, \
                                                    humanonlybam=self.humanOnlyBam, humanbetterbam=self.humanBetterBam, \
                                                    mouseonlybam=self.mouseOnlyBam, mousebetterbam=self.mouseBetterBam, \
                                                    humanlossbam=self.humanLossBam, mouselossbam=self.mouseLossBam), dependencies=[SortBamN_1, SortBamN_2])
        else:
            BamCmp = self.addTask("BamCmp", dependencies=[SortBamN_1, SortBamN_2])

        #step five
        if self.step < 5:
            SortBam = self.addTask("SortBam", self.bamsort(core=10, mem="3G", sortbamfile=self.humanBetterSortBam, bamfile=self.humanBetterBam,\
                                                       tmpdir=self.tmpdir), dependencies=[BamCmp])
        else:
            SortBam = self.addTask("SortBam",  dependencies=[BamCmp])

        # step six
        if self.step < 6:
            markDuplicate = self.addTask("markduplicate", self.Markduplicate(tmpdir=self.tmpbasedir, input_bam=self.humanBetterSortBam,
                                                                             marked_duplicated_bam=self.mkdupsortbam,
                                                                             marked_dup_metrics_txt=self.mkmatrix),
                                         dependencies=[SortBam])
        else:
            markDuplicate = self.addTask("markduplicate", dependencies=[SortBam])

        # step_seven
        if self.step < 7:
            BamIndex = self.addTask("BamIndex", self.sortbamindex(self.mkdupsortbam), dependencies=[markDuplicate])
        else:
            BamIndex = self.addTask("BamIndex", dependencies=[markDuplicate])

        # step_eight
        if self.step < 8:
            mutect2 = self.addTask("mutect2", self.Mutect2(tmpdir=self.tmpbasedir, ref=self.ref, inputbam=self.mkdupsortbam, rg=self.sample_id,
                                                           vcffile=self.vcffile), dependencies=[BamIndex])
        else:
            mutect2 = self.addTask("mutect2", dependencies=[BamIndex])


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument("-c", "--core", type=int, default=30)
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
    wf = bamcmp_workflow(fq1=args.fastq_1, fq2=args.fastq_2, tmpdir=args.tmp, ref=args.ref, workdir=args.wkdir, resultdir=args.output_dir, step=args.step)
    retval = wf.run(mode="local", nCores=60)
    sys.exit(retval)


