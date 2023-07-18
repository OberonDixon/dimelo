import dimelo as dm
import pandas as pd
import numpy as np
import sqlite3
import matplotlib.pyplot as plt
import pysam
import time
import argparse
import os

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

bam_filename = os.path.basename(bam_path)

start_time = time.time()
(pile_coordinates,valid_pile_dict,modified_pile_dict) = dm.parse_bam(
    bam_path,
    'test',output_dir,
    basemods=('N:A+m:N','N:C+m:G'),
    region=region,
    thresholds=[204,204],
    referenceGenome=genome_path)
end_time = time.time()
print(end_time-start_time)
np.savez(file=f'{output_dir}/{bam_filename}_{region}_piles.npz',
         pile_coordinates=pile_coordinates,
         valid_CpG=valid_pile_dict['N:C+m:G'],
        modified_CpG=modified_pile_dict['N:C+m:G'],
        valid_A=valid_pile_dict['N:A+m:N'],
        modified_A=modified_pile_dict['N:A+m:N'])