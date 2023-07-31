import pyBigWig
from dimelo.utils.file_saver import SubregionData,ProcesswiseTaskBuilder
import numpy as np

def save_subregion_batch(
    subregion_data: SubregionData,
    outDir_temp: str,
    format_key: str,
):
    temp_file_prefix = f'{subregion_data.string}_temp_{format_key}'
    
    chr_sizes = {
        subregion_data.chromosome:subregion_data.end
    }
    
    rd_bw = pyBigWig.open(f"{outDir_temp}/{temp_file_prefix}_read_depth.bw","w")
    
    rd_bw.addHeader([(k, v) for k, v in chr_sizes.items()])
    
    rd_bw.addEntries(
        np.repeat(subregion_data.chromosome,len(subregion_data.pile_coordinates)),
        subregion_data.pile_coordinates,
        subregion_data.pile_coordinates+1,
        subregion_data.read_depth
    )
    
    rd_bw.close()
    
    for basemod in subregion_data.basemods:
        mod_bw = pyBigWig.open(f"{outDir_temp}/{temp_file_prefix}_mod{basemod}.bw","w")        
        mod_bw.addHeader([(k, v) for k, v in chr_sizes.items()])        
        mod_bw.addEntries(
            np.repeat(subregion_data.chromosome,len(subregion_data.pile_coordinates)),
            subregion_data.pile_coordinates,
            subregion_data.pile_coordinates+1,   
            subregion_data.modified_pile_dict[basemod],
        )        
        mod_bw.close()
        
        valid_bw = pyBigWig.open(f"{outDir_temp}/{temp_file_prefix}_valid{basemod}.bw","w")
        valid_bw.addHeader([(k, v) for k, v in chr_sizes.items()])       
        valid_bw.addEntries(
            np.repeat(subregion_data.chromosome,len(subregion_data.pile_coordinates)),
            subregion_data.pile_coordinates,
            subregion_data.pile_coordinates+1,   
            subregion_data.valid_pile_dict[basemod],
        )       
        valid_bw.close()
        
        mod_over_valid_bw = pyBigWig.open(f"{outDir_temp}/{temp_file_prefix}_mod-over-valid{basemod}.bw","w")
        mod_over_valid_bw.addHeader([(k, v) for k, v in chr_sizes.items()])
        mod_over_valid_bw.addEntries(
            np.repeat(subregion_data.chromosome,len(subregion_data.pile_coordinates)),
            subregion_data.pile_coordinates,
            subregion_data.pile_coordinates+1,   
            np.where(subregion_data.valid_pile_dict[basemod] != 0,
                     subregion_data.modified_pile_dict[basemod]/
                     subregion_data.valid_pile_dict[basemod],
                     0),                      
        )        
        mod_over_valid_bw.close()
           
def merge_subregion_batches(
    processwise_tasks: ProcesswiseTaskBuilder,
    sampleName: str,
    outDir: str,
    outDir_temp: str,
    format_key: str,
    basemods: list,
):
    merged_file_path = f'{outDir}/{sampleName}_{format_key}_read_depth.bw'
    target_bw = pyBigWig.open(merged_file_path,"w")
    target_bw.addHeader([(k, v) for k, v in processwise_tasks.chr_ranges.items()])
        
    for single_process_tasks in processwise_tasks.core_assignments.values():
        for subregion_tasks_list in single_process_tasks.tasks.values():
            for subregion_task in subregion_tasks_list:
                temp_file_path = f'{outDir_temp}/{subregion_task.string}_temp_{format_key}_read_depth.bw'
                source_bw = pyBigWig.open(temp_file_path)
                chrom = subregion_task.chromosome
                
                starts,ends,values = map(list,zip(*source_bw.intervals(chrom)))
                
                target_bw.addEntries(
                    np.repeat(chrom,len(starts)),
                    starts,
                    ends,
                    values,
                )
    
    target_bw.close()
        
    for basemod in basemods:
        for track_type in ['mod','valid','mod-over-valid']:
            merged_file_path = f'{outDir}/{sampleName}_{format_key}_{track_type}{basemod}.bw'
            target_bw = pyBigWig.open(merged_file_path,"w")
            target_bw.addHeader([(k, v) for k, v in processwise_tasks.chr_ranges.items()])
            for single_process_tasks in processwise_tasks.core_assignments.values():
                for subregion_tasks_list in single_process_tasks.tasks.values():
                    for subregion_task in subregion_tasks_list:
                        temp_file_path = f'{outDir_temp}/{subregion_task.string}_temp_{format_key}_{track_type}{basemod}.bw'
                        source_bw = pyBigWig.open(temp_file_path)
                        chrom = subregion_task.chromosome
                        
                        starts,ends,values = map(list,zip(*source_bw.intervals(chrom)))
                        
                        target_bw.addEntries(
                            np.repeat(chrom,len(starts)),
                            starts,
                            ends,
                            values,
                        )

def create_bw(
    path:str,
    chr_sizes: dict,
    begin,
    end,
    values,
):
    pass