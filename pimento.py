
from pathlib import Path
import click

from bin.thresholds import MIN_STD_PRIMER_THRESHOLD
from bin.standard_primer_matching import get_primer_props, parse_std_primers, save_out

@click.group()
def cli():
    pass


@cli.command(
    "std",
    options_metavar="-i <fastq/fastq.gz> -p <primers_dir> -o <output_prefix>",
    short_help="Perform the standard primer strategy for primer inference"
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
    default=Path("./standard_primers")
)
@click.option(
    "-m",
    "--minimum_primer_threshold", 
    help="The minimum proportion of reads a standard primer has to be present in to be considered in inference. Default value of 0.60.", 
    type=float, 
    default=MIN_STD_PRIMER_THRESHOLD
)
@click.option(
    "-o", "--output_prefix", required=True, help="Prefix to output file.", type=str
)
def standard_primer_strategy(input_fastq: Path, primers_dir: Path, minimum_primer_threshold: float, output_prefix: str) -> None:

    std_primer_dict_regex, std_primer_dict = parse_std_primers(
        primers_dir
    )  # Parse std primer library into dictionaries
    results = get_primer_props(
        std_primer_dict_regex, input_fastq, minimum_primer_threshold
    )  # Find all the std primers in the input and select most common
    save_out(results, output_prefix, std_primer_dict)

if __name__ == "__main__":
    cli()