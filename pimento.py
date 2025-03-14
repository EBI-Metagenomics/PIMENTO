from pathlib import Path
import click

from bin.standard_primer_matching import (
    get_primer_props,
    parse_std_primers,
    write_std_output,
)
from bin.are_there_primers import atp_in_this_sample, write_atp_output
from bin.generate_bcv import generate_bcv_for_single_strand, write_bcv_output

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
    res_df.to_csv(f"{output_prefix}_mcp_cons.tsv", sep="\t")


if __name__ == "__main__":
    cli()
