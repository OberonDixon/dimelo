import pandas as pd
from typing import List, Union
import pysam
from collections import defaultdict
import numpy as np
import pandas as pd
import sys
from pympler import asizeof

class Region(object):
    def __init__(self, region: Union[str, pd.Series]):
        """Represents a region of genetic data.
        Attributes:
                - chromosome: string name of the chromosome to which the region applies
                - begin: integer start position of region
                - end: integer end position of region
                - size: length of region
                - string: string representation of region
                - strand: string specifying forward or reverse strand; either "+" or "-" (default +)
        """
        self.chromosome = None
        self.begin = None
        self.end = None
        self.size = None
        self.string = None
        self.strand = "+"

        if isinstance(region, str):  # ":" in region:
            # String of format "{CHROMOSOME}:{START}-{END}"
            try:
                self.chromosome, interval = region.replace(",", "").split(":")
                try:
                    # see if just integer chromosomes are used
                    self.chromosome = int(self.chromosome)
                except ValueError:
                    pass
                self.begin, self.end = [int(i) for i in interval.split("-")]
            except ValueError:
                raise TypeError(
                    "Invalid region string. Example of accepted format: 'chr5:150200605-150423790'"
                )
            self.size = self.end - self.begin
            self.string = f"{self.chromosome}_{self.begin}_{self.end}"
        elif isinstance(region, pd.Series):
            # Ordered sequence containing [CHROMOSOME, START, END] and optionally [STRAND], where STRAND can be either "+" or "-"
            self.chromosome = region[0]
            self.begin = region[1]
            self.end = region[2]
            self.size = self.end - self.begin
            self.string = f"{self.chromosome}_{self.begin}_{self.end}"
            # strand of motif to orient single molecules
            # if not passed just keep as all +
            if len(region) >= 4:
                if (region[3] == "+") or (region[3] == "-"):
                    self.strand = region[3]
                # handle case of bed file with additional field that isn't strand +/-
                else:
                    self.strand = "+"
            else:
                self.strand = "+"
        else:
            raise TypeError(
                "Unknown datatype passed for Region initialization"
            )
            
class SubregionTask(Region):
    VALID_POSITIONS = {'first', 'internal', 'last', 'only'}
    
    def __init__(
        self, 
        fileName:str,
        region: Union[str, pd.Series, Region], 
        intraregion_position: str,
        task_bases:int,
        formats_list: list,
    ):
        if isinstance(region,str) or isinstance(region,pd.Series):
            super().__init__(region)
        elif isinstance(region,Region):
            # If a Region object is passed in, copy its attributes
            self.chromosome = region.chromosome
            self.begin = region.begin
            self.end = region.end
            self.size = region.size
            self.string = region.string
            self.strand = region.strand
        
        if intraregion_position not in self.VALID_POSITIONS:
            raise ValueError(f"Invalid intraregion_position. Accepted values are {self.VALID_POSITIONS}")
        
        self.intraregion_position = intraregion_position
        self.bam_file = fileName
        self.task_bases = task_bases
        self.formats_list = formats_list

class SubregionData(SubregionTask):
    def __init__(
        self,
        subregion_task: SubregionTask,
        basemods: list,
    ):
        
        self.__dict__.update(vars(subregion_task))
        self.basemods = basemods
        
        if self.begin>=self.end:
            print('begin:end:',self.begin,self.end)
        
        self.read_depth = np.zeros(self.end-self.begin)
        
        self.modified_pile_dict = {}
        self.valid_pile_dict = {}
        self.pile_coordinates = np.arange(self.begin,self.end)
        
        for basemod_identifier in basemods:
            self.modified_pile_dict[basemod_identifier] = np.zeros(self.end-self.begin)
            self.valid_pile_dict[basemod_identifier] = np.zeros(self.end-self.begin)
        
        self.read_lists_dict = defaultdict(list)
            
    def add_read(
        self,
        read_name,
        read_chr,
        read_pos,
        read_is_forward,
        ref_seq,
        read_seq_aligned,
        reference_coordinates,
        read_coordinates,
        valid_coordinates_list_disambiguated,
        modified_coordinates_list_disambiguated,
    ):
        # Process individual read for reads dict
        if self.intraregion_position in ['internal','last'] and read_pos<self.begin:
            # This read has already been saved, skip it
            pass
        else:
            self.read_lists_dict['read_name'].append(read_name)
            self.read_lists_dict['chr'].append(read_chr)
            self.read_lists_dict['pos'].append(read_pos)
            self.read_lists_dict['is_forward'].append(read_is_forward)
            self.read_lists_dict['ref_seq'].append(ref_seq)
            self.read_lists_dict['read_seq_aligned'].append(read_seq_aligned)
            self.read_lists_dict['read_coordinates'].append(read_coordinates)
            self.read_lists_dict['reference_coordinates'].append(
                reference_coordinates)
            
            for basemod_index,basemod_identifier in enumerate(self.basemods):
                self.read_lists_dict[f'val{basemod_identifier}'].append(
                    valid_coordinates_list_disambiguated[basemod_index])
                self.read_lists_dict[f'mod{basemod_identifier}'].append(
                    modified_coordinates_list_disambiguated[basemod_index])
                
#             if self.end - read_pos < 1000:
#                 read_lists_dict_size = 0
#                 for key,value in self.read_lists_dict.items():
#                     key_size = asizeof.asizeof(value)
#                     if type(value[0])==np.ndarray:
#                         print(key,key_size,'type:',value[0].dtype,len(value),'elements, last one ',len(value[-1]),'slong')
#                     else:
#                         print(key,key_size,'type:',type(value[0]))
#                     read_lists_dict_size+=key_size

#                 print('reads memory usage:',read_lists_dict_size)
        
        # Perform pileup operation
        # create masks for valid reference_coordinates
        valid_mask = (reference_coordinates >= self.pile_coordinates[0]) & (reference_coordinates <= self.pile_coordinates[-1])
        read_indices = np.searchsorted(self.pile_coordinates,reference_coordinates[valid_mask])
        
        self.read_depth[read_indices]+=read_coordinates[valid_mask]
        for basemod_index,basemod_identifier in enumerate(self.basemods):
            self.modified_pile_dict[basemod_identifier][read_indices]+=(modified_coordinates_list_disambiguated
                                                                   [basemod_index][valid_mask])
            self.valid_pile_dict[basemod_identifier][read_indices]+=(valid_coordinates_list_disambiguated
                                                                [basemod_index][valid_mask])

class SingleProcessTasks:
    def __init__(
        self,
        process_id: int,
        formats_list: list,
    ):
        self.process_id = process_id
        self.tasks = defaultdict(list)
        self.total_bases = 0
    def add_task(
        self,
        region: Region,
        subregion_task: SubregionTask,
    ):
        self.tasks[region.string].append(subregion_task)
        self.total_bases += subregion_task.task_bases
class ProcesswiseTaskBuilder:
    def __init__(
        self, 
        region_list: List[Region], 
        fileName: str,
        formats_list: list,
        num_cores: int, 
        mem_allowance: int,
    ):
        self.region_list = region_list
        self.formats_list = formats_list
        self.num_cores = num_cores
        
        self.chr_ranges = {}
 
        self.construct_loadbalanced_tasks(
            fileName=fileName,
            batch_num_bases=mem_allowance/num_cores
        )
    
    def construct_loadbalanced_tasks(
        self, 
        fileName: str,
        batch_num_bases: int,
    ):     
        # Map out read across regions
        
        bam = pysam.AlignmentFile(fileName,"rb",check_sq=True)
        
        regionwise_read_positions = {}
        regionwise_read_lengths = {}
        
        for region in self.region_list:
            if region.chromosome not in self.chr_ranges:
                self.chr_ranges[region.chromosome]=region.end
            elif self.chr_ranges[region.chromosome]<region.end:
                self.chr_ranges[region.chromosome]=region.end
            positions = []
            readlens = []
            for read in bam.fetch(reference=region.chromosome,start=region.begin,end=region.end):
                positions.append(read.pos)
                readlens.append(read.query_alignment_length)
            regionwise_read_positions[region]=positions
            regionwise_read_lengths[region]=readlens
            
        total_read_bases = sum([sum(readlens) for readlens in regionwise_read_lengths.values()])
        approx_bases_per_core = (total_read_bases//self.num_cores)+1
        
        # Build single process task lists
        self.core_assignments = {}
        core_index = 0
        self.core_assignments[core_index] = SingleProcessTasks(core_index,self.formats_list)
                
        # We will add reads to a batch until we reach batch_num_bases, then split to a new subregion
        bases_in_batch = 0
        # We will add subregions to a core until we reach this value, then go to the next core
        bases_assigned_to_core = 0
        
        # Iterate through the positions and readlens for reads identified in each window
        for (region,positions,readlens) in [(region,regionwise_read_positions[region],regionwise_read_lengths[region]) 
                                            for region in self.region_list]:
            # The first subregion will have type_label 'first', meaning reads saved from these subregions 
            # are truncated at the start so they don't go out-of-coordinate-range
            subregion_start = region.begin
            type_label = 'first'
            # Iterate through the position and readlen values for reads in this window
            for position,readlen in zip(positions,readlens):
                # If we have exceeded the batch size, define a new subregion
                if bases_in_batch > batch_num_bases:
                    subregion_end = max(subregion_start+1,position-1)
                    if(subregion_start+1>position-1):
                        print(f'Warning: too little memory, batch {bases_in_batch} at {region.chromosome} hit limit {batch_num_bases} for core {core_index}')
                    self.core_assignments[core_index].add_task(
                        region,
                        SubregionTask(
                            fileName,
                            pd.Series([region.chromosome,subregion_start,subregion_end]),
                            type_label,
                            bases_in_batch,
                            self.formats_list,
                        )
                    )
                    subregion_start = position
                    print('bases in batch (batch size exceeded)',bases_in_batch)
                    bases_in_batch = 0
                # If we have exceeded the per-core approximate limit, define a new subregion and a new
                # SingleProcessTasks entry
                if bases_assigned_to_core > approx_bases_per_core and position!=last_position:
                    subregion_end = max(subregion_start,position-1)
                    if(subregion_start+1>position-1):
                        print(f'Warning: end of allotment mismatch, batch {bases_in_batch} at {region.chromosome} start {subregion_start} end {position-1}')
                    if bases_in_batch>0:
                        self.core_assignments[core_index].add_task(
                            region,
                            SubregionTask(
                                fileName,
                                pd.Series([region.chromosome,subregion_start,subregion_end]),
                                type_label,
                                bases_in_batch,
                                self.formats_list,
                            )
                        )
                    subregion_start = position
                    core_index += 1
                    self.core_assignments[core_index] = SingleProcessTasks(core_index,self.formats_list)
                    print('bases in batch (core allotment reached)',bases_in_batch)
                    bases_assigned_to_core = 0
                    bases_in_batch = 0
                    # All subregions after the first per window are 'internal', meaning reads
                    # they save should never be truncated on either end
                    type_label = 'internal'
                # Add the bases from the current read to the total for the core
                # This is after the check so that we err on the side of more bases per core,
                # ensuring we never assign more cores than we have (i.e. we will assign a bit too
                # much to all the cores except the last, which is ok)
                bases_assigned_to_core += readlen
                bases_in_batch += readlen
                last_position = position
            if subregion_start < region.end:
                if subregion_start == region.begin:
                    # If a window is not split at all into different subregions, its reads get
                    # truncated at both the beginning and the end
                    type_label = 'only'
                else:
                    # The last subregion of a region has reads it saves truncated at the end
                    type_label = 'last'
                if bases_in_batch>0:
                    self.core_assignments[core_index].add_task(
                        region,
                        SubregionTask(
                            fileName,
                            pd.Series([region.chromosome,subregion_start,region.end]),
                            type_label,
                            bases_in_batch,
                            self.formats_list,
                        )
                    )
                print(f'bases in batch ({type_label} for core)',bases_in_batch)
        
        