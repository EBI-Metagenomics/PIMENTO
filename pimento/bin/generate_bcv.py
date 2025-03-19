#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2025 EMBL - European Bioinformatics Institute
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

from collections import defaultdict
from pathlib import Path

import pandas as pd
import numpy as np

from pimento.bin.pimento_utils import (
    get_read_count,
    compute_windowed_base_conservation,
    build_list_of_base_counts,
    fetch_read_substrings,
)
from pimento.bin.thresholds import MAX_READ_COUNT, BCV_WINDOW_SIZE


def generate_bcv_for_single_strand(path: Path, rev: bool = False) -> defaultdict[float]:
    """
    Generate base-conservation vectors in a stepwise and windowed manner for a fastq file.

    For a continuous range of starting indices (2 to 25), generate substrings of window size of 5 bases.
    Calculate the average conservation of the most common base at each index of a window.
    The resulting list of base-conservations can be considered a conservation vector and used to
    identify inflection points where the conservation suddenly changes.

    Output a dictionary where:
        key -> an index starting point e.g. base 10
        val -> the avg conservation of the most common base for the substring window from base 10 to 15 (inclusive)
    """

    res_dict = defaultdict(float)
    start_range = range(2, 25, 1)  # Range of starting indices

    for start in start_range:

        read_count = get_read_count(
            path, file_type="fastq"
        )  # get read count for fastq file

        max_line_count = 0
        if read_count > MAX_READ_COUNT:
            max_line_count = MAX_READ_COUNT

        read_substring_count_dict = fetch_read_substrings(
            path, BCV_WINDOW_SIZE, rev, start, max_line_count
        )  # get substring count dict
        base_counts = build_list_of_base_counts(
            read_substring_count_dict, BCV_WINDOW_SIZE
        )  # list of base conservation dicts for substring
        base_conservation, cons_seq = compute_windowed_base_conservation(
            base_counts, read_count, max_line_count=max_line_count
        )  # get list of max base conservations for each index

        res_dict[start] = np.mean(base_conservation)  # compute the mean

    return res_dict


def write_bcv_output(
    fwd_out: defaultdict[float] = {}, rev_out: defaultdict[float] = {}
):
    """
    Generate Pandas dataframe out of bcv dictionary.

    Output looks like this (when both F and R are requested):
        2	3	4
    F	0.7814975041597337	0.8736772046589019	0.9434276206322796
    R	0.9010981697171381	0.9082861896838601	0.90369384359401

    Columns are the starting indices. Row labels are the strand.
    """

    final_res_dict = defaultdict(list)
    df_ind = []

    # Check if fwd strand was requested
    if fwd_out:
        [final_res_dict[key].append(fwd_out[key]) for key in fwd_out.keys()]
        df_ind.append("F")

    # Check if rev strand was requested
    if rev_out:
        [final_res_dict[key].append(rev_out[key]) for key in rev_out.keys()]
        df_ind.append("R")

    res_df = pd.DataFrame.from_dict(final_res_dict)
    res_df.index = df_ind

    return res_df
