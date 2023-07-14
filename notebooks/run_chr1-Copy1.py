import dimelo as dm
import pandas as pd
import numpy as np
import sqlite3
import matplotlib.pyplot as plt
import pysam
import time

combined_bam_filepath = '/clusterfs/nilah/oberon/downloads/prod_CTCF_winnowmap_guppy_merge.sorted.bam'
megalodon_bam_filepath = '/clusterfs/nilah/oberon/dimelo_dev_source/test_inputs/deep_ctcf_mod_mappings_merge.sorted.bam'
two_color_r9 = '/clusterfs/nilah/oberon/dimelo_dev_source/test_inputs/lmnb1-accessibility_20220214_megalodon/barcode10_rabbit-abcam/mod_mappings.10.sorted.bam'
test_r10 = '/clusterfs/nilah/oberon/dimelo_dev_source/test_inputs/r10/Dorado_R10_calls.bam'
hp1_bam_filepath = '/clusterfs/nilah/oberon/downloads/phased/ALLCTCF_guppy_winnowmap_merge_chr11_NanoMethPhase_HP1.bam'
hp2_bam_filepath = '/clusterfs/nilah/oberon/downloads/phased/ALLCTCF_guppy_winnowmap_merge_chr11_NanoMethPhase_HP2.bam'
pacbio_bam_filepath = '/clusterfs/nilah/oberon/dimelo_dev_source/test_inputs/PacBio/GM12878.CENPC.cell1.winnowmap.sorted.bam'
output_dir = '/clusterfs/nilah/oberon/dimelo_dev_source/test_outputs/'
sql_output = 'test_output'
hp1_sample_name = 'gm12878_ctcf_hp1'
genome_path = '/clusterfs/nilah/oberon/dimelo_dev_source/test_inputs/genomes/chm13.draft_v1.0.fasta'
v1_1_genome_path = '/clusterfs/nilah/oberon/genomes/chm13.draft_v1.1.fasta'
hg38_genome_path = '/clusterfs/nilah/ayesha/basenji/baselines/gm12878/data/hg38.fa'

region = 'chr1:0-250000000'
start_time = time.time()
(pile_coordinates,valid_pile_dict,modified_pile_dict) = dm.parse_bam(
    megalodon_bam_filepath,
    'test',output_dir,
    basemods=('N:A+m:N','N:C+m:G'),
    region=region,
    thresholds=[204,204],
    referenceGenome=genome_path)
end_time = time.time()
print(end_time-start_time)
np.savez(file=f'/clusterfs/nilah/oberon/dimelo_dev_source/test_outputs/{region}_piles_204.npz',
         pile_coordinates=pile_coordinates,
         valid_CpG=valid_pile_dict['N:C+m:G'],
        modified_CpG=modified_pile_dict['N:C+m:G'],
        valid_A=valid_pile_dict['N:A+m:N'],
        modified_A=modified_pile_dict['N:A+m:N'])