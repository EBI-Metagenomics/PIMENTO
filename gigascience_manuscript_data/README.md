## GigaScience Manuscript Data

This README describes the relationships between the data in this directory and the figures that were generated out of this data for the GigaScience manuscript and its supplementary material. The data is first divided into subdirectories for both the manuscript figures and the supplementary figures, with each containing subdirectories related to figures, containing the figures themselves and the data they were generated from.

```bash
.
├── manuscript_figures
│   ├── figures2-3-4-5
│   │  ├── are_there_primers_results.tar.gz
│   │  ├── marine_outcomes_plot.png
│   │  ├── marine_primer_counts.png
│   │  ├── marine_primer_proportions.png
│   │  ├── marine_primer_region_counts.png
│   │  └── std_primer_results.tar.gz
│   ├── figure6
│   │   ├── are_there_primers_results.tar.gz
│   │   ├── soil_outcomes_plot.png
│   │   └── std_primer_results.tar.gz
│   ├── figure7
│   │   ├── are_there_primers_results.tar.gz
│   │   ├── humangut_outcomes_plot.png
│   │   └── std_primer_results.tar.gz
│   └── figure8
│       ├── primer_identity_boxplots.png
│       └── std_auto_primer_comparison.csv
│  
└── supplementary_figures
    ├── S1
    │   ├── benchmarking_runtimes_std.csv
    │   ├── benchmarking_std_runtime_vs_file_size.png
    │   └── benchmarking_std_runtime_vs_read_count.png
    ├── S2
    │   ├── benchmarking_auto_runtime_vs_file_size.png
    │   ├── benchmarking_auto_runtime_vs_read_count.png
    │   └── benchmarking_runtimes_auto.csv
    └── S3
        ├── benchmarking_are_there_primers_runtime_vs_file_size.png
        ├── benchmarking_are_there_primers_runtime_vs_read_count.png
        └── benchmarking_runtimes_are_there_primers.csv
```

## manuscript_figures/
The manuscript has eight figures, seven of which were generated from the data in this directory (the final figure is a diagram not a plot). This directory has four subdirectories, one for each source of data making up at least one figure. Figures 2, 3, 4, and 5 are made up of the same PIMENTO outputs, and so share a subdirectory.

### figures2-3-4-5/
Figures 2 to 5 represent the PIMENTO analyses related to the Marine Survey Dataset. This subdirectory contains two compressed `.tar.gz` files that contain PIMENTO outputs for both the `are_there_primers` utility, and the standard primer matching method, with one output for each run in the dataset.

#### are_there_primers_results.tar.gz
When uncompressed, this archive contains files named like `ERR10176043_general_primer_out.txt`. The files contain two lines, each representing an orientation (5'-3' and 3'-5'), which can have a 0 or 1 as a value. A 1 represents a positive prediction for a primer existing at that orientation, with a negative prediction represented by a 0. These files were used to generate the `G` and `g` outcomes for Figure 2.

#### std_primer_results.tar.gz
When uncompressed, this archive contains files named like `ERR10176043_std_primer_out.txt`. The file is either empty if no standard primer was found, or contains up to three lines that look like this:
```
18S
1560F: 0.990624577742207
2035R: 0.9890135018553358
```
If at least one standard primer is found, the first line contains the predicted amplified region for that primer or primer pair (in this case `18S`). Following this line, this file contains one line for each primer found, with a maximum of one per orientation (in this case both `1560F` and `2035R`). These lines also contain the proportion of reads a primer was found in (in this case both `0.990624577742207` and `0.9890135018553358`). Figures 3, 4, and 5 were generated from the data inside these files.

### figure6
Figure 6 of the manuscript was generated in the exact same way as Figure 2, but for the Soil Survey Dataset, and therefore has the same two compressed archives.

### figure7
Figure 7 of the manuscript was generated in the exact same way as Figure 2, but for the Human Gut Survey Dataset, and therefore has the same two compressed archives.

### figure8
Figure 8 was generated from the single `.csv` file contained in this subdirectory, `std_auto_primer_comparison.csv`. This file looks like this:

```
pair,direction,std_file,auto_file,std_seq_id,auto_seq_id,std_len,auto_len,alignment_score,alignment_length,matches,mismatches,gaps,pct_identity_aln,pct_identity_vs_std,std_seq,auto_seq
ERR12749352,forward,ERR12749352_std_primers.fasta,ERR12749352_auto_primers.fasta,341F,F_auto,17,19,3.0,21,12,3,6,57.14,70.59,CCTACGGGNGGCWGCAG,DDDDCCTACGGGDKGCAGC
ERR12749352,reverse,ERR12749352_std_primers.fasta,ERR12749352_auto_primers.fasta,806BR,R_auto,20,18,6.0,23,14,1,8,60.87,70.0,GGACTACNVGGGTWTCTAAT,DDDDGACTACNVGGGTAT
ERR14753993,forward,ERR14753993_std_primers.fasta,ERR14753993_auto_primers.fasta,341F,F_auto,17,18,13.0,18,16,1,1,88.89,94.12,CCTACGGGNGGCWGCAG,CCTACGGGDGGCWGCAGT
ERR14753993,reverse,ERR14753993_std_primers.fasta,ERR14753993_auto_primers.fasta,805R,R_auto,21,21,19.0,21,20,1,0,95.24,95.24,GACTACHVGGGTATCTAATCC,GACTACHRGGGTATCTAATCC
```

Specifically, the boxplots of Figure 8 were generated from using the `direction`, `pct_identity_aln`, and `pct_identity_vs_std` columns, but the file contains the predicted primer sequences for extra visibility of the alignments.

## supplementary_figures
The supplementary material has three sections with figures (S1, S2, and S3), which are all about benchmarking the runtime performance of the three different PIMENTO utilities based on file size and read count. Each section has its own subdirectory.

### S1
The figures in S1 were generatedfrom the single `.csv` file contained in this subdirectory, `benchmarking_runtimes_std.csv`. This file looks like this:

```
file,real_time_sec,user_time_sec,sys_time_sec,file_size,read_count
SRR2657613_1.fastq.gz,95.37,293.98,3.86,1.51 GB,"13,316,864"
SRR2657613_2.fastq.gz,91.33,275.41,4.56,1.40 GB,"13,316,864"
ERR9597096_1.fastq.gz,38.07,110.70,1.79,532.84 MB,"8,947,245"
ERR9597096_2.fastq.gz,42.56,122.58,2.15,563.96 MB,"8,947,245"
```

Specifically, the bar plots were generated from the `real_time_sec`, `file_size`, and `read_count` columns.

### S2
The figures in S2 were generated in the exact same way as in S1, but using the `benchmarking_runtimes_auto.csv` file.

### S3
The figures in S3 were generated in the exact same way as in S1, but using the `benchmarking_runtimes_are_there_primers.csv` file.
