# Concatenated Influenza Workflow


## concat_seqs.py Script

This script is designed to concatenate influenza sequences together in a strict, pre-defined order on a per-file basis.
Alternatively, it can also extract an individual segment of your choosing per-file.  

### Input

The script receives either a single FASTA file or a folder of FASTA files as input. 

- If it receives a single file, it will produce a single file as output.
- If it receives a folder, it will produce a folder as output.

### Workflows

- `--concat`: will concatenate all 8 flu segments together on a per-file basis and ensure they are always in consistent order. The pre-defined order is their numbering scheme of #1-8 (PB2, PB1, PA, HA, NP, NA, M, NS). If fewer than 8 segments are present, the script will skip the input file.
- `--segment [SEGMENT]`: will extract a single flu segment that you specify. 


## Reference Folder

This folder contains a set of scripts used to curate a concatenated reference for use in alignments and trees. 
The workflow will create 1) a concatenated FASTA file and 2) a custom GFF3 file containing gene annotations for this sequence, all while ensuring consistency between both files. 