# PIMENTO
<img src="PIMENTO.png" width="500">

A **P**r**IME**r infere**N**ce **TO**olkit to facilitate large-scale calling of metabarcoding amplicon sequence variants.

## How PIMENTO works
PIMENTO’s employs a dual primer inference strategy, which are:

- **Standard primer search**: based on fuzzy regex search queries to a library of curated standard primer sequences.
- **Primer cutoff prediction**: based on the identification of the primer cutoff point from analysis of patterns of base-conservation at the beginning (and end, for single-end libraries) of reads. Consensus sequences are then generated as inferred primers using the predicted cutoff.

PIMENTO also implements an "are there primers?" function to predict the presence of primers in sequencing reads in case no standard primer was found. This method is helpful in cases where it isn't known whether primer sequences are still present in the reads, and checking manually would not be trivial, i.e. for large-scale analysis pipelines.

## How to install

PIMENTO is available on [PyPi](https://pypi.org/project/mi-pimento/). To install it from PyPi with `pip` just run:

`pip install mi-pimento`

PIMENTO is also available on bioconda and can be installed like this with conda/mamba:

`conda install -c bioconda mi-pimento`

## How to run

![PrimerInferenceWorkflow](PrimerInferenceWorkflow.png)

You can run either PIMENTO strategy with a single command. The tool will look for primers on either end, so both strategies will work on paired-end, single-end, or merged paired-end sequencing reads (though you would have to run it twice unmerged paired-end sequencing reads, one for each end).

```
pimento --help
Usage: pimento [OPTIONS] COMMAND [ARGS]...

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  are_there_primers     Predict whether primers are present in the input reads
  auto                  Perform the primer cutoff strategy for primer
                        inference
  choose_primer_cutoff  Choose the optimal primer cutoff point.
  find_cutoffs          Find potential cutoffs using a BCV output.
  gen_bcv               Generate the base-conservation vector(s) (BCV)
  std                   Perform the standard primer strategy for primer
                        inference
```


### Standard primer matching

To run the standard primer strategy:
```bash
pimento std -i <fastq/fastq.gz> -p <primers_dir> -o <output_prefix>
```

#### Inputs

`-i <fastq/fastq.gz>`: the input FASTQ reads file.

`-p <primers_dir>`: the path to the standard primers library to be used, with the default being PIMENTO's library. You can use your own library, or extend PIMENTO's. If using a different library than the default, make sure the primer FASTA files have this format:

```
>341F
CCTACGGGNGGCWGCAG
>338F
ACTCCTACGGGAGGCAGCA
>805R
GACTACHVGGGTATCTAATCC
>785R
CTACCAGGGTATCTAATCC
```

Where forward strand primers have the character `F` as the final character, and vice versa `R` for reverse strand primers.

`-o <output_prefix>`: the prefix to be used on output files.

#### Outputs

`<output_prefix>_std_primers.fasta`: FASTA file containing the best found single or pairs of primers. Empty if none were found.

`<output_prefix>_std_primer_out.txt`: Text file containing the read proportions of the best found primers.

`all_standard_primer_proportions.txt`: Text file logging all the read proportions for every single searched primer.

### Primer cutoff prediction

To run the primer cutoff strategy:
```bash
pimento auto -i <fastq/fastq.gz> -st [FR/F/R] -o <output_prefix>
```

NB: Running `pimento auto` executes the three subcommands `generate_bcv`, `find_cutoffs`, `choose_primer_cutoff` sequentially. You can therefore run each step of this workflow individually if you wish.

#### Inputs

`-i <fastq/fastq.gz>`: the input FASTQ reads file.

`-st [FR/F/R]`: the selection of strands to perform primer inference for - F for forward, R for reverse, FR for both.

`-o <output_prefix>`: the prefix to be used on output files.

#### Outputs

`<output_prefix>_auto_primers.fasta`: FASTA file containing the inferred primer sequences using the predicted optimal cutoffs.

### Are there primers?

To run the "are there primers?" utility:
```bash
pimento are_there_primers -i <fastq/fastq.gz> -o <output_prefix>
```

#### Inputs

`-i <fastq/fastq.gz>`: the input FASTQ reads file.

`-o <output_prefix>`: the prefix to be used on output files.


#### Outputs

`<output_prefix>_general_primer_out.txt`: Text file containing a 1 or 0 depending on if a primer was found on the forward strand (first line) and the reverse strand (second line).
