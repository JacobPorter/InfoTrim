#!/usr/bin/python
import InfoTrim
import datetime
import sys
import optparse
import subprocess
from Utilities import SeqIterator

"""
@author: Jacob Porter
TODO: test
nohup InfoTrim_automate.py -s -a -r 100 ./infotrim.automate.simulate.test.input ./infotrim.automate.simulate.test.output &> infotrim.automate.simulate.test.err &
"""

def automate(input_file, output_file_prefix, process_amount, accuracy, simulate, read_amount):
    process_amount = str(process_amount)
    if simulate:
        #subprocess.call(["~/DWGSIM/dwgsim", "-e", "0.012", "-E", "0.012", "-d", "250", "-s", "30", "-S", "0", "-N", str(read_amount), "-c", "2", "-1", "200", "-2", "200", "-f", "TACGTACGTCTGAGCATCGATCGATGTACAGC",  "/research/jsporter/Data/genome/GRCh38.p9/methyl-convert/GRCh38.p9.methyl-convert.GA.fa", input_file + ".GA"])
        subprocess.call(["/home/jsporter/DWGSIM/dwgsim", "-e", "0.012", "-E", "0.012", "-d", "250", "-s", "30", "-S", "0", "-N", str(read_amount), "-c", "0", "-1", "100", "-2", "100",  "/research/jsporter/Data/genome/GRCh38.p9/methyl-convert/GRCh38.p9.methyl-convert.GA.fa", input_file + ".GA"])
        subprocess.call(["/home/jsporter/DWGSIM/dwgsim", "-e", "0.012", "-E", "0.012", "-d", "250", "-s", "30", "-S", "0", "-N", str(read_amount), "-c", "0", "-1", "100", "-2", "100",  "/research/jsporter/Data/genome/GRCh38.p9/methyl-convert/GRCh38.p9.methyl-convert.CT.fa", input_file + ".CT"])
        subprocess.call(["/home/jsporter/BS_Simulation/dwgsim_postprocess", "-p", "2", "-o", "1", input_file + ".GA" + ".bwa.read2.fastq"], stdout = open(input_file + ".GA" + ".post", 'w'))
        subprocess.call(["/home/jsporter/BS_Simulation/dwgsim_postprocess", "-p", "1", "-o", "0", input_file + ".CT" + ".bwa.read1.fastq"], stdout = open(input_file + ".CT" + ".post", 'w'))
        input_file_2 = input_file + ".dwgsim.automate" + ".fastq"
        subprocess.call(["cat", input_file + ".CT" + ".post", input_file + ".GA" + ".post"], stdout = open(input_file_2, 'w'))
        input_file = input_file_2
    total_reads = SeqIterator.SeqIterator(input_file, file_type='fastq').count()
    subprocess.call(["/home/jsporter/walt-1.0/bin/walt", "-t", process_amount, "-i", "/research/jsporter/Data/genome/GRCh38.p9/walt/GRCh38.p9.multiLine.fa.index.walt.dbindex", "-o", input_file + ".walt.sam", "-r", input_file])
    subprocess.call(["/home/jsporter/bismark_v0.16.3/bismark", "-p", process_amount, "--sam",  "/research/jsporter/Data/genome/GRCh38.p9/bismark/", input_file])
    subprocess.call(["/home/jsporter/bwa-meth/bwameth.py", "--reference", "/research/jsporter/Data/genome/GRCh38.p9/bwa-meth/GRCh38.p9.fa", input_file], stdout = open(input_file + ".BWAMeth.sam", 'w'))
    subprocess.call(["countBWA.py", input_file + ".BWAMeth.sam"], stdout = open(input_file + ".BWAMeth.count", 'w'))
    bispin_args = ["BisPin_align.py", "-W", "-n", process_amount, "-i", "1", "-I", "2", "/research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa", input_file + ".BisPin.sam", input_file] 
    subprocess.call(bispin_args, stdout = open(input_file + ".BisPin.out", 'w'))
    if accuracy:
        subprocess.call(["calculateSimulationAccuracy.py", "-d", input_file + ".BisPin.sam", str(total_reads)], stdout = open(input_file + ".BisPin.acc", 'w'))
        subprocess.call(["calculateSimulationAccuracy.py", "-d", input_file.replace(".fastq", "") + "_bismark_bt2.sam", str(total_reads)], stdout = open(input_file + ".bismark.acc", 'w'))
        subprocess.call(["calculateSimulationAccuracy.py", "-d", input_file + ".BWAMeth.sam", str(total_reads)], stdout = open(input_file + ".BWAMeth.acc", 'w'))
        subprocess.call(["calculateSimulationAccuracy.py", "-d", input_file + ".walt.sam", str(total_reads)], stdout = open(input_file + ".walt.acc", 'w'))      
    for i in range(0, 11):
        for j in range(0, 11):
            if i + j > 10:
                break
            beta = i/10.0
            gamma = j/10.0
            alpha = (10 - i - j)/10.0
            m = 25
            t = 0
            #print alpha, beta, gamma
            output_file = "%s.s-%s.r-%s.m-%s.fastq" % (output_file_prefix, str(beta), str(gamma), str(m))
            InfoTrim.trim_reads(input_file, output_file, alpha, beta, gamma, m, t, process_amount=int(process_amount))
            subprocess.call(["/home/jsporter/bismark_v0.16.3/bismark", "-p", process_amount, "--sam",  "/research/jsporter/Data/genome/GRCh38.p9/bismark/", output_file])
            subprocess.call(["/home/jsporter/bwa-meth/bwameth.py", "--reference", "/research/jsporter/Data/genome/GRCh38.p9/bwa-meth/GRCh38.p9.fa", output_file], stdout = open(output_file + ".BWAMeth.sam", 'w'))
            subprocess.call(["countBWA.py", output_file + ".BWAMeth.sam"], stdout = open(output_file + ".BWAMeth.count", 'w'))
            subprocess.call(["/home/jsporter/walt-1.0/bin/walt", "-t", process_amount, "-i", "/research/jsporter/Data/genome/GRCh38.p9/walt/GRCh38.p9.multiLine.fa.index.walt.dbindex", "-o", output_file + ".walt.sam", "-r", output_file])
            bispin_args = ["BisPin_align.py", "-W", "-n", process_amount, "-i", "1", "-I", "2", "/research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa", output_file + ".BisPin.sam", output_file] 
            subprocess.call(bispin_args, stdout = open(output_file + ".BisPin.out", 'w'))
            if accuracy:
                subprocess.call(["calculateSimulationAccuracy.py", "-d", output_file + ".BisPin.sam", str(total_reads)], stdout = open(output_file + ".BisPin.acc", 'w'))
                subprocess.call(["calculateSimulationAccuracy.py", "-d", output_file.replace(".fastq", "") + "_bismark_bt2.sam", str(total_reads)], stdout = open(output_file + ".bismark.acc", 'w'))
                subprocess.call(["calculateSimulationAccuracy.py", "-d", output_file + ".BWAMeth.sam", str(total_reads)], stdout = open(output_file + ".BWAMeth.acc", 'w'))
                subprocess.call(["calculateSimulationAccuracy.py", "-d", output_file + ".walt.sam", str(total_reads)], stdout = open(output_file + ".walt.acc", 'w'))
        
def main():       
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <FASTQ reads file location> <output file prefix> "
    version = "0.1.0"
    p = optparse.OptionParser(usage = usage, version = version)
    p.add_option('--processes', '-p', help='Number of processes to use. [default: %default]', default=7)
    p.add_option('--read_amount', '-r', help='Approximate number of reads to simulate. [default: %default]', default=500000)
    p.add_option('--simulate', '-s', help="Simulate the reads.", action="store_true", default=False)
    p.add_option('--accuracy', '-a', help='Check the accuracy of DWGSIM reads', action="store_true", default=False)
    options, args = p.parse_args()
    if len(args) != 2:
        p.print_help()
        p.error("There must be two arugments.  Check the usage.")
    automate(args[0], args[1], int(options.processes), options.accuracy, options.simulate, int(options.read_amount))
    later = datetime.datetime.now()
    sys.stderr.write("Elapsed time for InfoTrim automate:\t%s\n" % (str(later - now)))
    
if __name__ == "__main__":
    main()
        