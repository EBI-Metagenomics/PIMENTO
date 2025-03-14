from pathlib import Path

import click
import pandas as pd

from bin.standard_primer_matching import (
    get_primer_props,
    parse_std_primers,
    write_std_output,
)
from bin.are_there_primers import atp_in_this_sample, write_atp_output
from bin.generate_bcv import generate_bcv_for_single_strand, write_bcv_output
from bin.find_cutoffs import find_bcv_inflection_points
from bin.choose_primer_cutoff import choose_cutoff_for_single_strand

from bin.thresholds import MIN_STD_PRIMER_THRESHOLD


@click.group()
def cli():
    pass


@cli.command(
    "std",
    options_metavar="-i <fastq/fastq.gz> -p <primers_dir> -o <output_prefix>",
    short_help="Perform the standard primer strategy for primer inference",
)
@click.option(
    "-i",
    "--input_fastq",
    required=True,
    help="Input fastq file to perform primer inference on.",
    type=click.Path(exists=True, path_type=Path, dir_okay=False),
)
@click.option(
    "-p",
    "--primers_dir",
    required=True,
    help="Input directory containing the standard primer library. Default uses the PIMENTO standard primer library.",
    type=click.Path(exists=True, path_type=Path, file_okay=False),
    default=Path("./standard_primers"),
)
@click.option(
    "-m",
    "--minimum_primer_threshold",
    help="The minimum proportion of reads a standard primer has to be present\
in to be considered in inference. Default value of 0.60.",
    type=float,
    default=MIN_STD_PRIMER_THRESHOLD,
)
@click.option(
    "-o", "--output_prefix", required=True, help="Prefix to output file.", type=str
)
def standard_primer_strategy(
    input_fastq: Path,
    primers_dir: Path,
    minimum_primer_threshold: float,
    output_prefix: str,
) -> None:

    std_primer_dict_regex, std_primer_dict = parse_std_primers(
        primers_dir
    )  # Parse std primer library into dictionaries
    results = get_primer_props(
        std_primer_dict_regex, input_fastq, minimum_primer_threshold
    )  # Find all the std primers in the input and select most common
    write_std_output(results, output_prefix, std_primer_dict)


@cli.command(
    "are_there_primers",
    options_metavar="-i <fastq/fastq.gz> -o <output_prefix>",
    short_help="Predict whether primers are present in the input reads",
)
@click.option(
    "-i",
    "--input_fastq",
    required=True,
    help="Input fastq file to predict the presence of primers for.",
    type=click.Path(exists=True, path_type=Path, dir_okay=False),
)
@click.option(
    "-o", "--output_prefix", required=True, help="Prefix to output file.", type=str
)
def are_there_primers(input_fastq: Path, output_prefix: str) -> None:

    fwd_primer_flag = atp_in_this_sample(
        input_fastq
    )  # Check for general primers in fwd
    rev_primer_flag = atp_in_this_sample(
        input_fastq, rev=True
    )  # Check for general primers in rev

    fwd_status = "0"
    rev_status = "0"
    # Flag for primer presence: 1 for yes 0 for no
    if fwd_primer_flag:
        print("Forward primer detected!")
        fwd_status = 1
    else:
        print("No forward primer detected")
    if rev_primer_flag:
        print("Reverse primer detected!")
        rev_status = 1
    else:
        print("No reverse primer detected")

    write_atp_output(
        (fwd_status, rev_status), output_prefix
    )  # Save primer flags to .txt file


@cli.command(
    "gen_bcv",
    options_metavar="-i <fastq/fastq.gz> -st [FR/F/R] -o <output_prefix>",
    short_help="Generate the base-conservation vector(s) (BCV)",
)
@click.option(
    "-i",
    "--input_fastq",
    required=True,
    help="Input fastq file to generate the BCV for.",
    type=click.Path(exists=True, path_type=Path, dir_okay=False),
)
@click.option(
    "-st",
    "--strand",
    help="The strand(s) to generate a BCV for. Values can be either F, R, or FR for both.",
    type=click.Choice(["FR", "F", "R"]),
    required=True,
)
@click.option(
    "-o", "--output_prefix", required=True, help="Prefix to output file.", type=str
)
def generate_base_conservation_vector(
    input_fastq: Path, strand: str, output_prefix: str
) -> None:

    res_df = ""

    # TODO: match-case statement is python 3.10>. We are currently locking the version
    # at version 3.9. The day we bump the version we should replace these if statements
    # with a match-case block.

    if strand == "FR":
        fwd_bcv = generate_bcv_for_single_strand(input_fastq)
        rev_bcv = generate_bcv_for_single_strand(input_fastq, rev=True)
        res_df = write_bcv_output(fwd_bcv, rev_bcv)
    elif strand == "F":
        fwd_bcv = generate_bcv_for_single_strand(input_fastq)
        res_df = write_bcv_output(fwd_bcv)
    elif strand == "R":
        rev_bcv = generate_bcv_for_single_strand(input_fastq, rev=True)
        res_df = write_bcv_output(rev_out=rev_bcv)

    # Save resulting dataframe to a tsv file
    res_df.to_csv(f"{output_prefix}_bcv.tsv", sep="\t")


@cli.command(
    "find_cutoffs",
    options_metavar="-i <bcv.tsv> -o <output_prefix>",
    short_help="Find potential cutoffs using a BCV output.",
)
@click.option(
    "-i",
    "--input_bcv",
    required=True,
    help="Input BCV file to identify potential cutoffs for.",
    type=click.Path(exists=True, path_type=Path, dir_okay=False),
)
@click.option(
    "-o", "--output_prefix", required=True, help="Prefix to output file.", type=str
)
def find_potential_cutoffs(input_bcv: Path, output_prefix: str):

    bcv_df = pd.read_csv(input_bcv, sep="\t", index_col=0)  # Read mcp_df
    inf_point_dict = find_bcv_inflection_points(
        bcv_df
    )  # Generate inflection points dict

    if len(inf_point_dict) > 0:  # If the inf_point_dict isn't empty..
        inf_point_df = pd.DataFrame.from_dict(
            inf_point_dict
        )  # .. turn it into a dataframe
        inf_point_df.to_csv(
            f"{output_prefix}_cutoffs.tsv", sep="\t", index=False
        )  # ..save it to a .tsv file

    else:  # If it is empty..
        fw = open(f"{output_prefix}_cutoffs.tsv", "w")  # ..make an empty file
        fw.close()


@cli.command(
    "choose_primer_cutoff",
    options_metavar="-i <fastq/fastq.gz> -p <cutoffs.tsv> -o <output_prefix>",
    short_help="Choose the optimal primer cutoff point.",
)
@click.option(
    "-i",
    "--input_fastq",
    required=True,
    help="Input fastq file to choose the optimal primer cutoff point for.",
    type=click.Path(exists=True, path_type=Path, dir_okay=False),
)
@click.option(
    "-p",
    "--primer_cutoffs",
    required=True,
    help="File containing the potential cutoff points to choose from.",
    type=click.Path(exists=True, path_type=Path, dir_okay=False),
)
@click.option(
    "-o", "--output_prefix", required=True, help="Prefix to output file.", type=str
)
def choose_primer_cutoff(input_fastq: Path, primer_cutoffs: Path, output_prefix: str):

    cutoffs_df = pd.read_csv(primer_cutoffs, sep="\t")

    f_slice = cutoffs_df[cutoffs_df.strand == "F"]  # get forward inflection points
    r_slice = cutoffs_df[cutoffs_df.strand == "R"]  # get reverse inflection points
    r_slice = r_slice.reset_index(drop=True)

    f_cutoff = ""
    r_cutoff = ""
    f_primer = ""
    r_primer = ""

    if not f_slice.empty:  # if there is a forward inflection point..
        cutoff_list = f_slice.inf_point.tolist()
        f_cutoff, f_primer = choose_cutoff_for_single_strand(
            input_fastq, cutoff_list
        )  # .. assess and select

    if not r_slice.empty:  # if there is a reverse inflection point..
        cutoff_list = r_slice.inf_point.tolist()
        r_cutoff, r_primer = choose_cutoff_for_single_strand(
            input_fastq, cutoff_list, rev=True
        )  # .. assess and select

    # Output cutoff point(s) to .txt file
    with open(f"{output_prefix}_chosen_cutoffs.txt", "w") as fw:
        if f_cutoff != "":
            fw.write(f"F: {f_cutoff}\n")
        if r_cutoff != "":
            fw.write(f"R: {r_cutoff}\n")

    # Output consensus primer sequence(s) to .fasta file
    with open(f"{output_prefix}_auto_primers.fasta", "w") as fw:
        if f_cutoff != "":
            fw.write(f">F_auto\n{f_primer}\n")
        if r_cutoff != "":
            fw.write(f">R_auto\n{r_primer}\n")


if __name__ == "__main__":
    cli()
