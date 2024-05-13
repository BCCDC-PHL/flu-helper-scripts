# Build Concatenated Reference

Use this set of scripts if you need to build a concatenated reference sequence. 
This process is more involved than simply concatenating sequences together because it assumes that you also require an associated GFF3 annotation file for your new reference.
The GFF3 annotation file is necessary for certain software programs to know the start/end positions of genomic regions within your concatenated reference sequence. 
Since no GFF3 file will exist for your new custom concatenated reference, these scripts will generate a custom one for you by fetching and the combining individual GFFs together. 

Here's how to use these scripts:

### Build Conda Environment

```
conda env create -f env.yml 
```

### Build MultiFASTA

Next, build a single multiFASTA file containing sequences in the order you want them to appear in the concatenated reference. 
Importantly, the headers in the multiFASTA file need to follow GenBank format, which is the accession number followed by a space before any other information (`>ACCESSION_NUMBER OTHER_INFO...`).
Here is an example for influenza A using a file called `flu_a_h1n1_reference.fasta`:

```
>NC_026438.1 Influenza A virus PB2
ATGCATCG
>NC_026435.1 Influenza A virus PB1
ATGCATCG
>NC_026437.1 Influenza A virus PA
ATGCATCG
>NC_026433.1 Influenza A virus HA
ATGCATCG
>NC_026436.1 Influenza A virus NP
ATGCATCG
>NC_026434.1 Influenza A virus NA
ATGCATCG
>NC_026431.1 Influenza A virus M
ATGCATCG
>NC_026432.1 Influenza A virus NS
ATGCATCG
```

### Fetch GFF3 Files From GenBank

This next step assumes the references you want to combine are accessible on GenBank.
In this case, this script minimizes the work needed to acquire the GFF3 files. 

```
conda activate concat
python fetch_gff3.py <MULTIFASTA> <OUTPUT_GFF_PATH>

# example
python fetch_gff3.py flu_a_h1n1_reference.fasta flu_gff_files
```

The output folder will contain GFF3 files with the filename format `<ACCESSION_NUMBER>.gff3`: 

```
flu_gff_files/
├─ NC_026438.1.gff3
├─ NC_026435.1.gff3
├─ NC_026437.1.gff3
├─ NC_026433.1.gff3
├─ NC_026436.1.gff3
├─ NC_026434.1.gff3
├─ NC_026431.1.gff3
├─ NC_026432.1.gff3
```
You can technically get your GFF3 files from any source.
The only requirement is that the name preceding `.gff3` matches the sequence name in the file. 

### Build Custom Reference Files

The final step is to generate the custom FASTA and GFF3 file. 
Using 1) your multiFASTA file and 2) the directory of GFF3 files you generated, run the second script:

```
python build_concat.py -f <MULTIFASTA> -g <OUTPUT_GFF_PATH> -o <CONCAT_GFF_FILENAME> -O <CONCAT_FASTA_FILENAME>

# example
python build_concat.py -f flu_a_h1n1_reference.fasta -g flu_gff_files -o flu_concat.gff3 -O flu_concat.fasta
```
