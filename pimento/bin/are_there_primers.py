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

from pathlib import Path
import numpy as np

from pimento.bin.pimento_utils import (
    get_read_count,
    compute_windowed_base_conservation,
    build_list_of_base_counts,
    fetch_read_substrings,
)

from pimento.bin.thresholds import ATP_WINDOW_SIZE, ATP_PREFIX_LENGTH


def atp_in_this_sample(input_fastq: Path, rev: bool = False) -> bool:
    """
    Predict the presence of primers based on windows of base conservation.

    Takes a fastq file as input. Extracts proportion of most common base for the first 100 bases.
    Computes the a threshold (Q3 - 0.15) based on this proportion and counts the number of bases below
    it in windows of 10 bases.
    If at least one of the first two windows contains at most one such a base, then the presence
    of a primer is flagged as true. A primer is also flagged as true if the combined count
    of bases below Q3 is at most 4.

    The output of this function is a boolean flag:
        True if a primer was identified
        False if a primer was not identified
    """

    read_count = get_read_count(
        input_fastq, file_type="fastq"
    )  # Get read count for fastq file

    read_substring_count_dict = fetch_read_substrings(
        input_fastq, ATP_PREFIX_LENGTH, rev=rev
    )  # substring dict where key is the substring and value is the count
    base_counts = build_list_of_base_counts(
        read_substring_count_dict, ATP_PREFIX_LENGTH
    )  # list of base conservation dicts for substrings
    base_conservation, cons_seq = compute_windowed_base_conservation(
        base_counts, read_count
    )  # get list of max base conservations for each index

    # Counter that will reset to 0 every 10 bases
    window_count = 0
    # Will append the window count to this list every 10 bases
    window_count_list = []
    # Compute Q3-based threshold
    max_cons = np.quantile(base_conservation, 0.75)
    threshold = max_cons - 0.15

    if max_cons < 0.75:
        threshold = 0.75
    # Immediately return false (no primer) if the max conservation is less than 0.6
    if max_cons < 0.6:
        return False

    # Loop through every base
    for counter, val in enumerate(base_conservation):
        if (
            counter % ATP_WINDOW_SIZE == 0 and counter != 0
        ):  # After looping through a window..
            window_count_list.append(window_count)  # ..append window count
            window_count = 0  # ..reset window count

        if (
            val < threshold
        ):  # If the conservation at i is less than threshold, increment count for the window
            window_count += 1

    primer_flag = False  # Initialise primer flag as false

    if (
        1 in window_count_list[:2] or 0 in window_count_list[:2]
    ):  # If window count is at most 1 of first two windows...
        primer_flag = True  # ..primer flag is true
    elif (
        sum(window_count_list[:2]) <= 4
    ):  # If sum of window counts of the first two windows is at most 4..
        primer_flag = True  # ..primer flag is true

    return primer_flag


def write_atp_output(results: tuple[bool, bool], output_prefix: str) -> None:
    """
    Save primer presence flags into output .txt file.

    1: primer exists
    0: primer doesn't exist

    First line will be the forward strand
    Second line will be the reverse strand
    """

    with open(f"{output_prefix}_general_primer_out.txt", "w") as fw:
        fw.write(f"{results[0]}\n")
        fw.write(f"{results[1]}\n")
