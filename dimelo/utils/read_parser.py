r"""
=================
readparser module
=================
.. currentmodule:: dimelo.readparser
.. autosummary::
    readpaser

readparser allows you to identify base modifications with context on .bam file reads, resolving ambiguities as needed

"""

import os
import json
import numpy as np
import pysam
from Bio.Seq import Seq
from typing import List, Tuple, Union

config_dicts = {}
config_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'basemods')
for filename in os.listdir(config_dir):
    if filename.endswith('.json'):
        file_path = os.path.join(config_dir,filename)

        with open(file_path,'r') as f:
            data = json.load(f)
            identifier = data.get('identifier',None)
            if identifier is not None:
                # Convert to a tuple of sets because that should be faster for running verifications later
                upstream_context_array_of_arrays = data['upstream_context']
                upstream_context_tuple_of_sets = tuple(set(sub_array) for sub_array in upstream_context_array_of_arrays)
                data['upstream_context'] = upstream_context_tuple_of_sets
                downstream_context_array_of_arrays = data['downstream_context']
                downstream_context_tuple_of_sets = tuple(set(sub_array) for sub_array in downstream_context_array_of_arrays)  
                data['downstream_context'] = downstream_context_tuple_of_sets
                config_dicts[identifier] = data                  
            else:
                print(f'Warning: basemod config {filename} is missing the "identifier" field')
CONFIGS = config_dicts

def parse_read_by_basemod(
    read:pysam.AlignedSegment,
    basemod_identifier:str,
    context_check_source:str='read',
    validate_with_reference:bool=False,
    genome:pysam.FastaFile=None,
    pipeline:str=None,
    threshold:int=128, 
) -> Tuple[np.ndarray,np.ndarray,np.ndarray,str]:
    """Pulls out arrays of modified and correctly contextualized bases for a given read and base modification

    Args:

    Returns:

     - reference coordinates array (length same as read length)
     - context-satisfying bases (length same as read length)
     - modified context-satisfying bases (length same as read length)
     - the string of the reference sequence

    """
    if basemod_identifier not in CONFIGS:
        print('error: basemod config is absent')

    config = CONFIGS[basemod_identifier]
    modified_bases = read.modified_bases
    basemod_key = None
    ########################################################################################
    # Identify basemod key
    ########################################################################################
    # The key format in the modified_bases dict from pysam is a tuple containing 
    # (modified_base,strand,modification_label)
    # This code currently does NOT check modification_label and should be modified to do so
    if modified_bases is not None:
        valid_key_counter = 0
        for key in modified_bases.keys():
#                 print(config['modified_base'])
#                 print(pipeline)
            if (key[config['modified_base_label_position']]==config['modified_base'] and
                (pipeline is None 
                or key[config['modification_label_position_by_pipeline'][pipeline]]
                    ==config['modification_label_by_pipeline'][pipeline])):
                basemod_key = key
#                     print(basemod_key)
                valid_key_counter += 1
        if valid_key_counter>1:
            print('error: more than one base modification label in .bam file. please specify pipeline and ensure .json config files contain appropriate label details.')
#             if valid_key_counter==0:
#                 print(list(modified_bases.keys()))
    ########################################################################################
    # Extract indices for modified bases that meet the threshold to be counted
    ########################################################################################
    # These indices tell us where, relative to the start of the read, we can find bases
    # that have been modified
    forward_sequence = read.get_forward_sequence()
    if forward_sequence is not None:
        valid_indices = range(0,len(read.get_forward_sequence()))
        read_indices = range(0,len(read.get_forward_sequence()))
        if basemod_key is not None:               
            if threshold>0:
                modified_indices = [coord_prob_tuple[0] for coord_prob_tuple 
                                    in modified_bases[basemod_key] if coord_prob_tuple[1]>=threshold]
            else:
                modified_indices = [coord_prob_tuple[0] for coord_prob_tuple 
                                    in modified_bases[basemod_key]]
        else:
            modified_indices = []
    else:
        valid_indices = []
        modified_indices = []
        read_indices = []
    ########################################################################################
    # Load up and process the read sequence for later comparisons, if needed
    ########################################################################################
    if context_check_source in ('read','both'):
        if read.is_forward:
            read_seq = read.query_sequence
        else:
            read_seq = str(Seq(read.query_sequence).complement())
#             print('read:',read_seq[0:10])
    else:
        read_seq = None


#         print('is_forward:',read.is_forward)



    ########################################################################################
    # Extract reference genome coordinates for read bases
    ########################################################################################
    # We always need this, because in order to have a single coordinate system for our pileups
    # we can't use the read coordinate system: reads have insertions and deletions sometimes
    reference_positions = read.get_reference_positions(full_length=True)
    read_start = min(coord for coord in reference_positions if coord is not None)
    read_end = max(coord for coord in reference_positions if coord is not None)+1
    reference_positions_rel = [position-read_start if position is not None else None for position in reference_positions]
    # If an index for a modified base has no corresponding reference genome coordinate,
    # we have no way of making that meaningful
    # The assumption here has to be that if the user wants their reads to get piled up or
    # peak called in a non-standard genome, they provide that .fasta for the upstream processing
    # steps i.e. the .bam file is aligned to it, and thus read.get_reference_positions() does not
    # return None for insertions
    valid_indices = [index for index in valid_indices if reference_positions[index] is not None]
    modified_indices = [index for index in modified_indices if reference_positions[index] is not None]
    read_indices = [index for index in read_indices if reference_positions[index] is not None]

    ########################################################################################
    # Load up the reference genome segment if needed
    ########################################################################################
    # We want to complement if on the reverse strand but NOT reverse complement, because our
    # reference_positions coordinates are not reversed but our context bases need to be complemented
    if validate_with_reference or context_check_source in ('reference','both'):
        fetched_seq = genome.fetch(read.reference_name,read_start,read_end)
        if read.is_forward:
            ref_seq = str(fetched_seq)
        else:                       
            ref_seq = str(Seq(str(fetched_seq)).complement())
    else:
        ref_seq = None

#         if ref_seq is not None:
#             print('ref:',ref_seq[0:10])

    ########################################################################################
    # Adjust indices to remove basemods that whose basecall doesn't match the reference
    ########################################################################################
    # If validate_with_reference is True, the index of bases that don't line up with the reference
    # is not counted as "modified" or "unmodified" it is simply removed from the list entirely
    modified_unvalidated_count = len(modified_indices)
    all_unvalidated_count = len(valid_indices)

    if validate_with_reference:
#             for index in valid_indices[1:-2]:
#                 print(basemod,ref_seq(reference_positions_rel[index-1:index+1]))
        modified_indices = [index for index in modified_indices 
                            if reference_positions_rel[index]<len(ref_seq) 
                            and ref_seq[reference_positions_rel[index]]
                                .upper()==config['modified_base']]
        valid_indices = [index for index in valid_indices 
                         if reference_positions_rel[index]<len(ref_seq) 
                         and ref_seq[reference_positions_rel[index]]
                             .upper()==config['modified_base']]

    modified_validated_count = len(modified_indices)
    all_validated_count = len(valid_indices)
#         print('validate?',validate_with_reference)
#         if validate_with_reference:
#             print('fraction that check out with reference',modified_validated_count/modified_unvalidated_count)

    ########################################################################################
    # Check upstream context
    ########################################################################################
    # If a base does not meet context requirements, its index is removed from the list
    # Context can be checked in the read itself, in the reference, or in both
    # We check one at a time which may not be optimal for multi-base context but is more readable
    # Upstream is reversed by convention because we are stepping our context distance upwards but
    # the upstream_context specifiers are in forward strand order
    for context_index,bases in enumerate(reversed(config['upstream_context'])):
        # If any base is valid, no more compute need be wasted
        if bases.issuperset(set(['A','T','C','G'])):
            continue
        else:
            if read.is_forward:
                distance = -context_index-1
            else:
                distance = context_index+1
            if context_check_source in ('reference','both'):
                modified_indices = [index for index in modified_indices 
                                    if len(reference_positions)>index+distance>0 
                                    and reference_positions[index+distance] is not None 
                                    and ref_seq[reference_positions_rel[index+distance]].upper() in bases]
                valid_indices = [index for index in valid_indices 
                                 if len(reference_positions)>index+distance>0 
                                 and reference_positions[index+distance] is not None 
                                 and ref_seq[reference_positions_rel[index+distance]].upper() in bases]
            if context_check_source in ('read','both'):
                modified_indices = [index for index in modified_indices 
                                    if len(read_seq)>index+distance>0 
                                    and read_seq[index+distance] in bases]
                valid_indices = [index for index in valid_indices 
                                 if len(read_seq)>index+distance>0 
                                 and read_seq[index+distance] in bases]


    ########################################################################################
    # Check downstream context
    ########################################################################################
    # If a base does not meet context requirements, its index is removed from the list
    # Context can be checked in the read itself, in the reference, or in both
    # We check one at a time which may not be optimal for multi-base context but is more readable               
    for context_index,bases in enumerate(config['downstream_context']):
        # If any base is valid, no more compute need be wasted
        if bases.issuperset(set(['A','T','C','G'])):
            continue
        else:
            if read.is_forward:
                distance = context_index+1
            else:
                distance = -context_index-1
            if context_check_source in ('reference','both'):
                modified_indices = [index for index in modified_indices 
                                    if 0<index+distance<len(reference_positions) 
                                    and reference_positions[index+distance] is not None 
                                    and ref_seq[reference_positions_rel[index+distance]].upper() in bases]
                valid_indices = [index for index in valid_indices 
                                 if 0<index+distance<len(reference_positions) 
                                 and reference_positions[index+distance] is not None 
                                 and ref_seq[reference_positions_rel[index+distance]].upper() in bases]
            if context_check_source in ('read','both'):
                modified_indices = [index for index in modified_indices 
                                    if 0<index+distance<len(read_seq) 
                                    and read_seq[index+distance] in bases]  
                valid_indices = [index for index in valid_indices 
                                 if 0<index+distance<len(read_seq) 
                                 and read_seq[index+distance] in bases]  

    modified_contextualized_count = len(modified_indices)
    all_contextualized_count = len(valid_indices)

    ########################################################################################
    # Build output arrays
    ######################################################################################## 

    reference_coordinates = np.arange(read_start,read_end,dtype=int)
    
    # Positions where the read aligns, so we can get read depth
    read_coordinates = np.zeros(read_end-read_start,dtype=int)
    if len(read_indices)>0:
        read_coordinates[np.array(reference_positions_rel)[np.array(read_indices)].astype(int)] = 1
    
    # Positions where the context is valid for a base modification, so we can get pileup denominator
    valid_coordinates = np.zeros(read_end-read_start,dtype=int)
    if len(valid_indices)>0:
        valid_coordinates[np.array(reference_positions_rel)[np.array(valid_indices)].astype(int)] = 1
    
    # Positions where the base-in-context was actually modified, so we can get pileup numerator. 
    # Report the actual prob if threshold was zero
    modified_coordinates = np.zeros(read_end-read_start,dtype=float)
    if len(modified_indices)>0:
        if threshold>0:
            modified_coordinates[np.array(reference_positions_rel)[np.array(modified_indices)].astype(int)] = 1
        else:
            modified_indices_set = set(modified_indices)
            filtered_modified_bases_tuples = [t for t in modified_bases[basemod_key] if t[0] in modified_indices_set]
            filtered_indices,filtered_values = zip(*filtered_modified_bases_tuples)
            modified_coordinates[np.array(reference_positions_rel)
                                 [np.array(filtered_indices)].astype(int)] = filtered_values/255
    
    # Index the read sequences to report back aligned with the valid and modified coordinates
    read_seq_aligned = "".join([read_seq[index] for index in range(len(read_seq)) 
                                if reference_positions[index] is not None])
    
    return (reference_coordinates,
            read_coordinates,
            valid_coordinates,
            modified_coordinates,
            ref_seq,
            read_seq_aligned,
           )

def resolve_basemod_ambiguities(
    basemods,
    valid_list,
    modified_list,    
) -> Tuple[list,list]:
    """Uses the .json config properties to guide decisions on treating base modifications with multiple valid contexts
    "
    " This version of resolve_basemod_ambiguities actually doesn't properly handle the fact that we *can* distinguish, 
    " in principle, between different modifications on the same base in the same context. It would be easy to fix that, 
    " the requisite information is already in the .json files, but with no use case or test case I decided not to implement
    " a solution at this time.
    """
    # Create a set of which possibilites there are for modified bases, i.e. which are the possible central
    # bases for varying contexts against which we want to check
    all_modified_base_options = set([config['modified_base'] for config in CONFIGS.values()])
    # Stack the valid coordinates list and modified coordinates list into numpy arrays for fast indexing
    valid_stack = np.stack(valid_list)
    modified_stack = np.stack(modified_list)

    # For each possible modified base, find which of the specified basemods are accessing that base and what coordinates
    # in the stacks corresponds to those basemods
    for modified_base in all_modified_base_options:
        stack_coordinates = np.array([coordinate for (coordinate,basemod) 
                             in enumerate(basemods) 
                             if CONFIGS[basemod]['modified_base']==modified_base],
                                     dtype=int)
        stack_coordinates_do_not_keep = np.array([coordinate for coordinate in stack_coordinates
                                                 if CONFIGS[basemods[coordinate]]['keep_if_ambiguous'] is False])

        # We only need to overwrite anything if at least one of the basemods centered on modified_base has keep_if_ambiguous
        # set to False, otherwise we can just skip to the next modified_base option
        if len(stack_coordinates_do_not_keep)>0:           
            sums = np.sum(valid_stack[stack_coordinates,:],axis=0)
            # This is technically not perfectly efficient if the reason sums>1 if not due to any of the
            # basemods within stack_coordinates_do_not_keep, but with only two basemod contexts per modifiable
            # base that will not occur, and this will still work regardless
            valid_stack[stack_coordinates_do_not_keep,sums>1] = 0
            modified_stack[stack_coordinates_do_not_keep,sums>1] = 0

    valid_list = list(valid_stack)
    modified_list = list(modified_stack)

    return (valid_list,modified_list)