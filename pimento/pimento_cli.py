from pathlib import Path
import importlib.metadata
from importlib.resources import files

import click
import pandas as pd
import pyfastx
from rich.console import Console
from rich import print

from pimento.bin.standard_primer_matching import (
    get_primer_props,
    parse_std_primers,
    write_std_output,
)
from pimento.bin.are_there_primers import atp_in_this_sample, write_atp_output
from pimento.bin.generate_bcv import generate_bcv_for_single_strand, write_bcv_output
from pimento.bin.find_cutoffs import find_bcv_inflection_points
from pimento.bin.choose_primer_cutoff import choose_cutoff_for_single_strand

from pimento.bin.thresholds import MIN_STD_PRIMER_THRESHOLD

console = Console()
__version__ = importlib.metadata.version("mi-pimento")
DEFAULT_STD_PRIMERS_PATH = files("pimento.standard_primers").joinpath("")


@click.group()
@click.version_option(__version__)
def cli():
    console.log(
        "[bold grey74]Starting new [bold green]PI[/bold green][bold red]MENTO[/bold red] run![bold grey74]\n"
    )
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
    default=Path(DEFAULT_STD_PRIMERS_PATH),
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
@click.option(
    "--merged",
    default=False,
    is_flag=True,
    help="Flag for running the standard primer strategy in `merged` mode, \
which is necessary if the input sequence file is made up of merged paired-end \
or single-end reads. Specifically, it will reverse the substrings that are searched\
for reverse primers, and use the complement of said reverse primers.",
)
def standard_primer_strategy(
    input_fastq: Path,
    primers_dir: Path,
    minimum_primer_threshold: float,
    output_prefix: str,
    merged: bool,
) -> None:
    """Runs the standard primer matching strategy for primer inference.
    A library of standard primers will be searched in input reads using fuzzy regex to identify
    primers with high degrees of proportions. The best single, or pair, of primers found will be
    outputted into a FASTA file (provided any were indeed found).
    PIMENTO comes with its own standard primer library which will be used as a default. Users can
    however give their own custom input library, or extend the default library to suit their needs.

    :param input_fastq: The input FASTQ file.
    :type input_fastq: Path
    :param primers_dir: Directory containing FASTA files of primer sequences for the library.
    The only format requirement is that forward strand primer names in the FASTA headers
    finish with the character `F`, and vice versa `R` for reverse strand primers. See PIMENTO's
    standard library for examples.
    :type primers_dir: Path
    :param minimum_primer_threshold: The minimum matching proportion threshold for a primer to be considered.
    The current default value for this threshold is a proportion of 0.60 of reads. Users can customise this value.
    :type minimum_primer_threshold: float
    :param output_prefix: The prefix to be used on output files.
    :type output_prefix: str
    """
    print(
        "[bold grey74]Running [bold green]standard primer strategy[/bold green].[/bold grey74]"
    )
    print(
        f"[bold grey74]Input FASTQ file: [bold green]{input_fastq}[/bold green][/bold grey74]"
    )
    print(
        f"[bold grey74]Standard primer library: [bold green]{primers_dir}[/bold green][/bold grey74]"
    )
    print(
        f"[bold grey74]Output prefix: [bold green]{output_prefix}[/bold green][/bold grey74]"
    )
    print("")

    with console.status("[bold yellow]Loading standard primer library..."):
        std_primer_dict_regex, std_primer_dict, primer_count = parse_std_primers(
            primers_dir, merged
        )  # Parse std primer library into dictionaries
        console.log(
            "[bold green]Loading standard primer library :white_check_mark:[/bold green]\n"
        )

    print(
        f"[bold grey74][bold green]{primer_count} primers[/bold green] loaded[/bold grey74]"
    )

    with console.status("[bold yellow]Searching for standard primers..."):
        results = get_primer_props(
            std_primer_dict_regex, input_fastq, minimum_primer_threshold, merged
        )  # Find all the std primers in the input and select most common
        console.log(
            "[bold green]Searching for standard primers :white_check_mark:[/bold green]\n"
        )

    std_primers_fasta, std_primers_info = write_std_output(
        results, output_prefix, std_primer_dict, merged
    )

    if results:
        primer_fasta = pyfastx.Fasta(str(std_primers_fasta), build_index=False)
        print("[bold grey74]Found standard library primers:[/bold grey74]")
        for primer in primer_fasta:
            print(f"[bold green]>{primer[0]}[/bold green]")
            print(f"[bold red]{primer[1]}[/bold red]")
        print("")
        print(
            f"[bold grey74]Final primers found saved to [bold green]{std_primers_fasta}[/bold green][/bold grey74]"
        )
        print(
            f"[bold grey74]Proportion information saved to [bold green]{std_primers_info}[/bold green][/bold grey74]"
        )
        print(
            "[bold grey74]All proportions saved to [bold green]all_standard_primer_proportions.txt[/bold green]"
            "[/bold grey74]"
        )
    else:
        print("No standard library primers found in reads!")
        print(
            "[bold grey74]Log of proportions for all primers found in "
            "[bold green]all_standard_primer_proportions.txt[/bold green][/bold grey74]"
        )


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
    """Checks for the presence of primers in input reads.
    Using patterns of base-conservation values for input reads, this function makes a prediction
    for whether there are primers at either ends of reads. The flag for whether primers are detected
    (1 for yes, 0 for no) for both strands is saved to an output file.

    :param input_fastq: The input FASTQ file.
    :type input_fastq: Path
    :param output_prefix: The prefix to be used on output files.
    :type output_prefix: str
    """
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
) -> Path:
    """Generates base-conservation vector(s) for input fastq files.
    To be used by the primer cutoff inference strategy, this function computes the base-conservation
    vector for input reads. Users can choose to perform this computation for either forward or reverse
    strands of the reads, or both.

    :param input_fastq: The input FASTQ file.
    :type input_fastq: Path
    :param strand: The strand(s) to perform primer inference for. Values can be either F, R, or FR for both.
    :type strand: str
    :param output_prefix: The prefix to be used on output files.
    :type output_prefix: str
    :return: TSV file containing the base-conservation vector.
    :rtype: Path
    """

    with console.status("[bold yellow]Generating base-conservation vector..."):

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

    console.log(
        "[bold green]Generating base-conservation vector :white_check_mark:[/bold green]"
    )

    return f"{output_prefix}_bcv.tsv"


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
def find_potential_cutoffs(input_bcv: Path, output_prefix: str) -> Path:
    """Find the potential cutoffs in a base-conservation vector output file.
    By computing negative inflection points in the base-conservation vectors,
    this function identifies potential cutoff points to be considered as primer cutoff
    points in the next and final step of the primer cutoff inference strategy. These potential
    cutoffs are saved to a TSV file.

    :param input_bcv: TSV file containing the base-conservation vector.
    :type input_bcv: Path
    :param output_prefix: The prefix to be used on output files.
    :type output_prefix: str
    :return: TSV file containing the potential cutoff points to consider.
    :rtype: Path
    """
    with console.status("[bold yellow]Finding potential cutoffs..."):

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

    console.log("[bold green]Finding potential cutoffs :white_check_mark:[/bold green]")

    return Path(f"{output_prefix}_cutoffs.tsv")


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
def choose_primer_cutoff(
    input_fastq: Path, primer_cutoffs: Path, output_prefix: str
) -> Path:
    """Choose the optimal primer cutoff point from an input set.
    Using patterns of base-conservation in reads and a set of potential cutoff points
    to choose from, this function makes a prediction for the "optimal" cutoff point
    to be used for primer inference. This function then generates a consensus sequence
    using the reads and the chosen cutoff point to be used as the final primer sequence.
    Performed for both strands, the final primer sequences are output into a FASTA file.

    :param input_fastq: The input FASTQ file.
    :type input_fastq: Path
    :param primer_cutoffs: TSV file containing the potential cutoff points to consider.
    :type primer_cutoffs: Path
    :param output_prefix: The prefix to be used on output files.
    :type output_prefix: str
    :return: The output FASTA file containing the inferred primer sequences.
    :rtype: Path
    """
    with console.status("[bold yellow]Choosing optimal primer cutoff point..."):

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

    console.log(
        "[bold green]Choosing optimal primer cutoff point :white_check_mark:[/bold green]\n"
    )

    return Path(f"{output_prefix}_auto_primers.fasta")


@cli.command(
    "auto",
    options_metavar="-i <fastq/fastq.gz> -st [FR/F/R] -o <output_prefix>",
    short_help="Perform the primer cutoff strategy for primer inference",
)
@click.option(
    "-i",
    "--input_fastq",
    required=True,
    help="Input fastq file to perform primer inference on.",
    type=click.Path(exists=True, path_type=Path, dir_okay=False),
)
@click.option(
    "-st",
    "--strand",
    help="The strand(s) to perform primer inference for. Values can be either F, R, or FR for both.",
    type=click.Choice(["FR", "F", "R"]),
    required=True,
)
@click.option(
    "-o", "--output_prefix", required=True, help="Prefix to output file.", type=str
)
@click.pass_context
def primer_cutoff_strategy(
    ctx: click.Context,
    input_fastq: Path,
    strand: str,
    output_prefix: str,
) -> None:
    """Runs the primer cutoff strategy for primer inference.
    Using patterns of base-conservation values for input reads, this strategy identifies potential cutoff
    points where a primer sequence could end. An optimal choice of cutoff based on the change in base-conservation
    before and after a cutoff is then made. Finally a consensus sequence is generated using this cutoff point
    and the input reads, which are saved to an output FASTA file.

    :param ctx: Click context - needed to invoke the multiple Click functions this "pipeline" uses.
    :type ctx: click.Context
    :param input_fastq: The input FASTQ file.
    :type input_fastq: Path
    :param strand: The strand(s) to perform primer inference for. Values can be either F, R, or FR for both.
    :type strand: str
    :param output_prefix: The prefix to be used on output files.
    :type output_prefix: str
    """

    print(
        "[bold grey74]Running [bold green]primer cutoff strategy[/bold green].[/bold grey74]"
    )
    print(
        f"[bold grey74]Input FASTQ file: [bold green]{input_fastq}[/bold green][/bold grey74]"
    )
    print(
        f"[bold grey74]Strands to perform strategy on: [bold green]{strand}[/bold green][/bold grey74]"
    )
    print(
        f"[bold grey74]Output prefix: [bold green]{output_prefix}[/bold green][/bold grey74]"
    )
    print("")

    bcv_df = ctx.invoke(
        generate_base_conservation_vector,
        input_fastq=input_fastq,
        strand=strand,
        output_prefix=output_prefix,
    )

    cutoffs_path = ctx.invoke(
        find_potential_cutoffs, input_bcv=bcv_df, output_prefix=output_prefix
    )

    inferred_primer_seq_path = ctx.invoke(
        choose_primer_cutoff,
        input_fastq=input_fastq,
        primer_cutoffs=cutoffs_path,
        output_prefix=output_prefix,
    )

    print("[bold grey74]Primer cutoff strategy finished![/bold grey74]")
    primer_fasta = pyfastx.Fasta(str(inferred_primer_seq_path), build_index=False)
    print("[bold grey74]Inferred primers:[/bold grey74]")
    for primer in primer_fasta:
        print(f"[bold green]>{primer[0]}[/bold green]")
        print(f"[bold red]{primer[1]}[/bold red]")
    print("")
    print(
        f"[bold grey74]Final primers found saved to [bold green]{inferred_primer_seq_path}[/bold green][/bold grey74]"
    )


if __name__ == "__main__":
    cli()
