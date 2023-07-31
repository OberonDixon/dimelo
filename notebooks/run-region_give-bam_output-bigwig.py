import dimelo as dm
import pandas as pd
import numpy as np
import sqlite3
import matplotlib.pyplot as plt
import pysam
import time
import argparse
import os
from pathlib import Path

parser = argparse.ArgumentParser(description="This script runs the dimelo parse_bam function on a whole chromosome with a specified bam input filepath.")

parser.add_argument("-b","--bam_path",help="path to load bam",type=str)
parser.add_argument("-r","--region",help="region to run",type=str)
parser.add_argument("-g","--genome",help="reference genome",type=str)
parser.add_argument("-o","--output_dir",help="output npz directory",type=str)

args = parser.parse_args()

bam_path = args.bam_path
region = args.region
genome_path = args.genome
output_dir = args.output_dir

bam_filename = Path(os.path.basename(bam_path)).stem

start_time = time.time()
dm.parse_bam(
    bam_path,
    bam_filename,output_dir,
    basemods=('mA','CpG'),
    region=region,
    thresholds=[204,204],
    referenceGenome=genome_path,
    checkAgainstRef=True,
    formats_list = ['bigwig'],
    cores=1,
    memory=100000000,
)
end_time = time.time()
print(end_time-start_time)