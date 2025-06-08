import pandas as pd
from Bio.Seq import reverse_complement
from .primer_design import add_upstream_sequence, add_downstream_sequence

class FlexProbeDict(pd.DataFrame):
    """
    A DataFrame subclass for handling flexible probe dictionaries.
    It allows for easy manipulation and access to probe data.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    
    def add_reverse_complement(self, sequence_column='target_sequence', rc_column='target_sequence_rc'):
        """Add reverse complement sequences to the DataFrame."""
        self[rc_column] = self[sequence_column].apply(reverse_complement)


#primer_design.add_upstream_sequence(
#    probe_dict,
#    capture_seq_series['Sequence_rc'],
#    input_column='target_sequence',
#    output_column=f"{capture_seq_series['Name']}_target_sequence",
#
#)


def add_upstream_capture_sequence(
    probe_dict: dict,
    capture_seq_name: str,
    capture_seq: str,
    input_column: str = 'target_sequence',
    output_column: str = 'capture_target_sequence',
    capture_seq_column: str = 'capture_sequence',
):
    '''Add a capture sequence for the flex probe library'''
    print(f"capture_seq_name: {capture_seq_name}")
    print(f"capture_seq: {capture_seq}")
    add_upstream_sequence(
        probe_dict,
        capture_seq,
        input_column=input_column,
        output_column=output_column,
    )
    # also add this capture sequence to the DataFrame
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            probe_dict[gk][tk][capture_seq_column] = capture_seq # save sequences
            probe_dict[gk][tk][f"{capture_seq_column}_name"] = capture_seq_name # also save names
            
def add_downstream_sequencing_primer_sequence(
    probe_dict: dict,
    primer_seq_name: str,
    primer_seq: str,
    input_column: str = 'capture_target_sequence',
    output_column: str = 'capture_target_sequence_sequencing_primer',
    primer_seq_column: str = 'sequencing_primer_sequence',
):
    '''Add a primer sequence for the flex probe library'''
    print(f"sequencing primer_seq_name: {primer_seq_name}")
    print(f"sequencing primer_seq: {primer_seq}")
    add_downstream_sequence(
        probe_dict,
        primer_seq,
        input_column=input_column,
        output_column=output_column,
    )
    # also add this primer sequence to the DataFrame
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            probe_dict[gk][tk][primer_seq_column] = primer_seq # save sequences
            probe_dict[gk][tk][f"{primer_seq_column}_name"] = primer_seq_name # also save names

def _full_flex_probe_basename(
    probe: pd.Series,
    capture_seq_column: str = 'capture_sequence',
    primer_seq_column: str = 'sequencing_primer_sequence',
):
    probe_name = [
        probe['gene_id'],
        probe['transcript_id'],
        f"shift:{probe['shift']}" ,
        f"GC:{probe['target_GC_left25']:.1f}-{probe['target_GC_right25']:.1f}",
        f"specificity:{probe['target_specificity']:.2f}" if 'target_specificity' in probe else '',
        f"isospecificity:{probe['target_isospecificity']:.2f}" if 'target_isospecificity' in probe else '',
        f"capture:{probe[capture_seq_column+'_name'].replace('_','-')}" if capture_seq_column+'_name' in probe else '',
        f"seqprimer:{probe[primer_seq_column+'_name'].replace('_','-')}" if primer_seq_column+'_name' in probe else '',
    ]
    return "_".join(probe_name)   

def _splitted_fled_probe_basename(
    probe: pd.Series,
    capture_seq_column: str = 'capture_sequence',
    primer_seq_column: str = 'sequencing_primer_sequence',
):
    if capture_seq_column in probe.keys():
        # this is left probe:
        probe_name = [
            probe['gene_id'],
            probe['transcript_id'],
            f"shift:{probe['shift']}" ,
            f"GC:{probe['target_GC_left25']:.1f}",
            f"specificity:{probe['target_specificity']:.2f}" if 'target_specificity' in probe else '',
            f"isospecificity:{probe['target_isospecificity']:.2f}" if 'target_isospecificity' in probe else '',
            f"capture:{probe[capture_seq_column+'_name'].replace('_','-')}" if capture_seq_column+'_name' in probe else '',
        ]
        return "_".join(probe_name) + "_left"
    elif primer_seq_column in probe.keys():
        # this is right probe:
        probe_name = [
            probe['gene_id'],
            probe['transcript_id'],
            f"shift:{probe['shift']}" ,
            f"GC:{probe['target_GC_right25']:.1f}",
            f"specificity:{probe['target_specificity']:.2f}" if 'target_specificity' in probe else '',
            f"isospecificity:{probe['target_isospecificity']:.2f}" if 'target_isospecificity' in probe else '',
            f"seqprimer:{probe[primer_seq_column+'_name'].replace('_','-')}" if primer_seq_column+'_name' in probe else '',
        ]
        return "_".join(probe_name) + "_right"
    else:
        raise ValueError("Probe does not contain capture or sequencing primer sequence columns, cannot determine probe name.")

def _split_flex_probe(
    probe: pd.Series,
    capture_seq_column: str = 'capture_sequence',
    primer_seq_column: str = 'sequencing_primer_sequence',
    target_seq_column: str = 'target_sequence',
    input_column: str = 'capture_target_sequence_sequencing_primer',
    split_pos: int = 25, # full probe length is 50, so split at 25
    output_key: str = 'Flex_probe_sequence',
    verbose: bool = False,
):
    # calculate lengths for left and right parts:
    left_length = len(probe[capture_seq_column]) + split_pos
    right_length = len(probe[primer_seq_column]) + (len(probe[target_seq_column]) - split_pos)
    if verbose:
        print(f"Splitting probe {probe.name} into left ({left_length}) and right ({right_length}) parts")
    # create left part
    left_probe = pd.Series({_k:_v for _k,_v in probe.items()
        if target_seq_column not in _k and primer_seq_column not in _k})
    left_probe[target_seq_column] = probe[target_seq_column][:split_pos]
    left_probe['original_'+target_seq_column] = probe[target_seq_column]
    left_probe[output_key] = probe[input_column][:left_length]
    # create right part:
    right_probe = pd.Series({_k:_v for _k,_v in probe.items()
        if target_seq_column not in _k and capture_seq_column not in _k})
    right_probe[target_seq_column] = probe[target_seq_column][split_pos:]
    right_probe['original_'+target_seq_column] = probe[target_seq_column]
    right_probe[output_key] = probe[input_column][-right_length:]
    # add name:
    
    # return both parts
    return left_probe, right_probe


def name_full_flex_probe_dict(
    probe_dict: dict, # flex probe dictionary to name, already assembled with capture and sequencing-primer sequences
    capture_seq_column: str = 'capture_sequence',
    primer_seq_column: str = 'sequencing_primer_sequence',
    input_column: str = 'capture_target_sequence_sequencing_primer',
    output_key: str = 'flex_probe_basename', # the key to use for the full probe name in the DataFrame
):
    """
    Name a flex probe dictionary with full probe names.
    
    Args:
        probe_dict (dict): The flex probe dictionary to name.
        
    Returns:
        dict: A dictionary of DataFrames, each keyed by the full probe name.
    """
    
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            if input_column not in probe_dict[gk][tk].columns:
                raise ValueError(f"Input column '{input_column}' not found in probe dictionary.")
            # name the full flex probe
            probe_dict[gk][tk][output_key] = pd.Series([_full_flex_probe_basename(
                _p,
                capture_seq_column=capture_seq_column,
                primer_seq_column=primer_seq_column,
                ) for _i, _p in probe_dict[gk][tk].iterrows()],
            )

def split_flex_probe_dict(
    probe_dict: dict, # flex probe dictionary to split, already assembled with capture and sequencing-primer sequences
    input_column: str = 'capture_target_sequence_sequencing_primer',
    capture_seq_column: str = 'capture_sequence',
    primer_seq_column: str = 'sequencing_primer_sequence',
    split_pos: int = 25, # full probe length is 50, so split at 25
    output_key: str = 'Flex_probe_sequence',
    ):
    """
    Split a flex probe dictionary into multiple DataFrames based on a specified column.
    
    Args:
        probe_dict (dict): The flex probe dictionary to split.
        split_by (str): The column name to split the DataFrame by.
        
    Returns:
        dict: A dictionary of DataFrames, each keyed by the unique values in the specified column.
    """

    
    left_probe_dict = {}
    right_probe_dict = {}
    
    for gk in probe_dict.keys():
        print(f"Processing Gene: {gk}")
        if gk not in left_probe_dict:
            left_probe_dict[gk] = {}
        if gk not in right_probe_dict:
            right_probe_dict[gk] = {}
        for tk in probe_dict[gk].keys():
            print(f"Processing Transcript: {tk}")
            if input_column not in probe_dict[gk][tk].columns:
                raise ValueError(f"Input column '{input_column}' not found in probe dictionary.")
            left_probes, right_probes = [],[]
            for _idx, _probe in probe_dict[gk][tk].iterrows():
                # split the probe
                left_probe, right_probe = _split_flex_probe(
                    _probe,
                    input_column=input_column,
                    split_pos=split_pos,
                    output_key=output_key,
                )
                # add name:
                left_probe[output_key+'_name'] = _splitted_fled_probe_basename(
                    left_probe,
                    capture_seq_column=capture_seq_column,
                    primer_seq_column=primer_seq_column,
                )
                right_probe[output_key+'_name'] = _splitted_fled_probe_basename(
                    right_probe,
                    capture_seq_column=capture_seq_column,
                    primer_seq_column=primer_seq_column,
                )
                # add to left and right lists
                left_probes.append(left_probe)
                right_probes.append(right_probe)
            # append:
            left_probe_dict[gk][tk] = pd.DataFrame(left_probes, index=probe_dict[gk][tk].index)
            right_probe_dict[gk][tk] = pd.DataFrame(right_probes, index=probe_dict[gk][tk].index)

    return left_probe_dict, right_probe_dict

def _select_probe_by_distance(
    probes: pd.DataFrame,
    distance: int = 0, # bases, 0 means right next to each other
    target_sequence_column: str = 'target_sequence',
):
    if 'shift' not in probes.columns:
        raise ValueError("Probes DataFrame must contain 'shift' column to select by distance.")
    
    shift_distance = distance + len(probes[target_sequence_column].tolist()[0])
    
    # version1: do the simplest:
    selected_probes = []
    current_shift = -1 * shift_distance # start with -1 to allow first probe to be selected
    for _, probe in probes.iterrows():
        if probe['shift'] >= current_shift + shift_distance:
            selected_probes.append(probe)
            current_shift = probe['shift']
    # convert to DataFrame
    selected_probes_df = pd.DataFrame(selected_probes)
    
    return selected_probes_df

def select_probe_dict_by_distance(
    probe_dict: dict,
    distance: int = 0, # bases, 0 means right next to each other
    target_sequence_column: str = 'target_sequence',
):
    """
    Select probes from a probe dictionary by distance.
    
    Args:
        probe_dict (dict): The probe dictionary to select from.
        distance (int): The distance in bases to select probes by.
        
    Returns:
        dict: A dictionary of DataFrames with selected probes.
    """
    sel_pb_dict = {}
    for _g in probe_dict:
        sel_pb_dict[_g] = {}
        for _t in probe_dict[_g]:
            sel_pb_dict[_g][_t] = _select_probe_by_distance(
                probe_dict[_g][_t],
                distance=distance,
                target_sequence_column=target_sequence_column,
            )
    
    return sel_pb_dict

def select_probe_dict_by_value(
    probe_dict: dict,
    value: str, # the value to select by, for example, isospecificity or specificity
    num_keep: int = 2, # number of probes to keep, if more than one probe has the same value, keep the first num_keep probes
    ascending: bool = True, # whether to keep low values (True) or high values (False)
):
    """
    Select probes from a probe dictionary by a specific value in a column.
    Args:
        probe_dict (dict): The probe dictionary to select from.
        value (str): The column name to select by.
        num_keep (int): Number of probes to keep for each gene and transcript.
        keep_low (bool): Whether to keep low values (True) or high values (False).
    Returns:
        dict: A dictionary of DataFrames with selected probes.
    """
    sel_pb_dict = {}
    for _g in probe_dict:
        sel_pb_dict[_g] = {}
        for _t in probe_dict[_g]:
            # sort by value
            sorted_probes = probe_dict[_g][_t].sort_values(by=value, ascending=ascending)
            # select top num_keep probes
            selected_probes = sorted_probes.head(num_keep)
            sel_pb_dict[_g][_t] = selected_probes.reset_index(drop=True)
    
    return sel_pb_dict
    