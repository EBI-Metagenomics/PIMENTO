#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2024 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
from collections import defaultdict
import os

from Bio.Seq import Seq
import regex

from mgnify_pipelines_toolkit.analysis.amplicon.amplicon_utils import (
    primer_regex_query_builder,
    get_read_count,
    fetch_mcp,
)


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=str,
        help="Path to merged FASTQ to look for primers",
    )
    parser.add_argument(
        "-p",
        "--primers",
        required=True,
        type=str,
        help="Path to directory containing standard primers fasta files",
    )
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output path")
    args = parser.parse_args()

    input = args.input
    primers = args.primers
    sample = args.sample
    output = args.output

    return input, primers, sample, output


def parse_std_primers(primers):
    """
    Parse the library of standard primers.

    Reads the fasta files in the given directory
    Primer names (which are the fasta headers) are labeled with F or R for 5'-3' and 3'-5' primers respectively

    Returns two dictionaries:
        std_primer_dict_regex
            key: region+primer name
            val: primer sequence from 5' to 3'
        std_primer_dict
            key: region+primer name
            val: primer sequence from 5' to 3' for forward primers, 3' to 5' for reverse
    """

    std_primer_dict_regex = defaultdict(defaultdict)
    std_primer_dict = defaultdict(defaultdict)

    dir = os.listdir(primers)
    dir = [f"{primers}/{path}" for path in dir]

    rev_flag = False

    for path in dir:
        region = path.split("/")[-1].split(".")[0]
        with open(path, "r") as fr:
            key = ""
            for line in fr:
                line = line.strip()
                if line[0] == ">":
                    if "R" in line:  # If a primer is a reverse primer
                        rev_flag = True
                    key = line[1:]
                else:
                    if rev_flag:
                        rev_conv = str(Seq(line).complement())
                        line = rev_conv
                        rev_flag = False

                    primer = primer_regex_query_builder(line)
                    std_primer_dict_regex[region][key] = primer
                    std_primer_dict[region][key] = line

    return std_primer_dict_regex, std_primer_dict


def run_primer_matching_once(input_path, input_primer, rev=False):
    """
    Run primer matching using the regex package.

    Takes one primer, strand, and fastq input
    Uses fuzzy matching to allow for at most one error (for sequencing errors)
    Returns number of reads matching given primer
    """

    match_count = 0.0

    mcp_count_dict = fetch_mcp(input_path, 50, rev=rev)

    for mcp in mcp_count_dict.keys():
        mcp = mcp.strip()
        res = regex.match(input_primer, mcp)
        if res is not None:
            match_count += mcp_count_dict[mcp]

    return match_count


def get_primer_props(std_primer_dict_regex, input_path):
    """
    Look for the standard primers in the input fastq file.

    Will loop through the dictionary of primers, using fuzzy regex matching to find matching primers.
    If a std primer is present above a set threshold proportion, it is collected. Both strands are searched for.
    If there is an std primer for both the F and R strands, the maximum prop for each strand is chosen and the pair
    is output as a combination.

    Returns a list containing two elements:
        max_region: the amplified region the chosen primers belong to
        max_primers: dictionary containing the F and/or R primers that were chosen
    """

    threshold = 0.60  # Arbitrary threshold for collecting a matched primer
    read_count = get_read_count(
        input_path, file_type="fastq"
    )  # Get read count of fastq file to calculate proportion with
    res_dict = defaultdict(defaultdict)

    # Loop through every primer region
    for region, primer in std_primer_dict_regex.items():
        res_dict[region]["F"] = {}
        res_dict[region]["R"] = {}

        # Loop through every primer of a certain region
        for primer_name, primer_seq in primer.items():

            region_name_str = f"{region};{primer_name}"
            primer_count = 0.0

            if "F" in primer_name:
                primer_count = run_primer_matching_once(
                    input_path, primer_seq, rev=False
                )  # Get proportion of a F primer with fuzzy regex matching
            elif "R" in primer_name:
                primer_count = run_primer_matching_once(
                    input_path, primer_seq, rev=True
                )  # Get proportion of a R primer with fuzzy regex matching

            try:
                primer_prop = primer_count / read_count
            except ZeroDivisionError:
                primer_prop = 0

            if "F" in primer_name:
                if (
                    primer_prop > threshold
                ):  # Only collect primer if it's above threshold
                    res_dict[region]["F"][primer_name] = primer_prop
            elif "R" in primer_name:
                if (
                    primer_prop > threshold
                ):  # Only collect primer if it's above threshold
                    res_dict[region]["R"][primer_name] = primer_prop

            print(f"{region_name_str}: {primer_prop}")

        # If an F or/and R primer wasn't found then just remove it from the dictionary
        if res_dict[region]["F"] == {}:
            res_dict[region].pop("F")
        if res_dict[region]["R"] == {}:
            res_dict[region].pop("R")

    singles = defaultdict(str)
    doubles = defaultdict(list)

    double_status = False  # Flag for whether primers were found on both strands

    #  Loop through every collected primer and put primers in singles or doubles
    for region in res_dict.keys():
        strands = res_dict[region]

        for strand in strands.keys():
            primers = strands[strand]
            max_prop = 0
            max_name = ""
            for primer_name, prop in primers.items():
                if prop > max_prop:
                    max_prop = prop
                    max_name = primer_name

            if len(strands.keys()) == 2:
                double_status = True
                doubles[region].append({max_name: max_prop})
            elif len(strands.keys()) == 1:
                singles[region] = {max_name: max_prop}

    max_region = ""
    max_primers = {}
    max_mean_prop = 0

    # if at least one pair of primers was collected
    if double_status:
        for (
            region
        ) in doubles:  # Loop through all pairs of primers and choose the best one
            primers = doubles[region]

            f_primer_name = list(primers[0].keys())[0]
            r_primer_name = list(primers[1].keys())[0]
            f_primer_prop = primers[0][f_primer_name]
            r_primer_prop = primers[1][r_primer_name]

            mean_prop = (f_primer_prop + r_primer_prop) / 2.0
            if mean_prop > max_mean_prop:
                max_mean_prop = mean_prop
                max_region = region
                max_primers = [
                    {f_primer_name: f_primer_prop},
                    {r_primer_name: r_primer_prop},
                ]

    else:
        for region in singles:  # Choose the best single primer
            primer = singles[region]
            primer_name = list(primer.keys())[0]
            prop = primer[primer_name]
            if prop > max_mean_prop:
                max_mean_prop = prop
                max_region = region
                max_primers = {primer_name: prop}

    if max_region == "":
        print("No standard library primers!")
        return []
    elif double_status:
        print("Standard library primers found!")
        print(f"Region: {max_region}")
        print(f"Forward Primer: {max_primers[0]}")
        print(f"Reverse Primer: {max_primers[1]}")

        return [max_region, max_primers[0], max_primers[1]]
    else:
        print("Standard library primer found on one strand!")
        print(f"Region: {max_region}")
        print(f"Primer: {max_primers}")

        return [max_region, max_primers]


def save_out(results, sample_id, output, std_primer_dict):
    """
    Save found std primers into a fasta file.
    """

    with (
        open(f"{output}/{sample_id}_std_primer_out.txt", "w") as fw_out,
        open(f"{output}/{sample_id}_std_primers.fasta", "w") as fw_seq,
    ):
        if results == []:
            fw_out.write("")
            fw_seq.write("")

        elif len(results) == 2:
            region = results[0]
            primer_name = list(results[1].keys())[0]
            primer_prop = results[1][list(results[1].keys())[0]]
            seq = std_primer_dict[region][primer_name]
            if "R" in primer_name:
                seq = str(Seq(seq).complement())
            fw_out.write(f"{region}\n")
            fw_out.write(f"{primer_name}: {primer_prop}")

            fw_seq.write(f">{primer_name}\n{seq}")

        elif len(results) == 3:
            region = results[0]
            f_primer_name = list(results[1].keys())[0]
            f_primer_prop = results[1][list(results[1].keys())[0]]
            f_seq = std_primer_dict[region][f_primer_name]
            r_primer_name = list(results[2].keys())[0]
            r_primer_prop = results[2][list(results[2].keys())[0]]
            r_seq = std_primer_dict[region][r_primer_name]
            r_seq = str(Seq(r_seq).complement())

            fw_out.write(f"{region}\n")
            fw_out.write(f"{f_primer_name}: {f_primer_prop}\n")
            fw_out.write(f"{r_primer_name}: {r_primer_prop}")

            fw_seq.write(f">{f_primer_name}\n{f_seq}\n")
            fw_seq.write(f">{r_primer_name}\n{r_seq}\n")


def main():

    input, primers, sample, output = parse_args()
    std_primer_dict_regex, std_primer_dict = parse_std_primers(
        primers
    )  # Parse std primer library into dictionaries
    results = get_primer_props(
        std_primer_dict_regex, input
    )  # Find all the std primers in the input and select most common
    save_out(results, sample, output, std_primer_dict)


if __name__ == "__main__":
    main()
