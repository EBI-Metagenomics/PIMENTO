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

from Bio.Seq import Seq
import numpy as np
import pandas as pd

from bin.amplicon_utils import (
    get_read_count,
    build_cons_seq,
    build_read_substring_cons_dict_list,
    fetch_read_substrings,
)
from bin.thresholds import MCP_MAX_LINE_COUNT


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=str,
        help="Path to fastq file to choose inflection point",
    )
    parser.add_argument(
        "-p", "--points", required=True, type=str, help="Path to inflection points file"
    )
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output path")

    args = parser.parse_args()

    path = args.input
    points = args.points
    sample = args.sample
    output = args.output

    return path, points, sample, output


def assess_inflection_point_mcp_for_sample(path, inf_point_list, rev=False):
    """
    Assess inflection point list, selecting one for automatic primer trimming.

    Takes as input a fastq file and a list of inflection points generated by "find_mcp_inflection_points.py".
    Computes the average conservation of mcp before inflection point and after.
    Gets the difference in avg. conservation between the pre- and post- points.
    Selects the inflection point with the maximum difference as the cutoff.
    If an inf point has a similar difference and is earlier than the max, we make a 'conservative' choice and
    replace it with the earlier cutoff.

    Returns the cutoff point and the consensus sequence 'forming' the automatically predicted primer
    """

    # TODO error handle for empty inflection point list

    start_confs = []  # pre-inf point conservations
    end_confs = []  # post-inf point conservations
    start_cons_lens = []  # list for storing lengths of pre-inflection point sequences
    cons_seq_list = []  # list for storing consensus sequences pre-inflection points

    do_not_include_list = [
        i + 5 for i in inf_point_list
    ]  # ignore conservation of inflection point in calculation

    read_count = get_read_count(path)  # get readcount from fastq

    max_line_count = 0
    if read_count > MCP_MAX_LINE_COUNT:
        max_line_count = MCP_MAX_LINE_COUNT

    n_prop = 0.8

    for start in inf_point_list:  # Looping through the pre-inflection point mcps
        mcp_len = start + 4  # length of pre-inf mcps is inflection point + 4

        mcp_count_dict = fetch_read_substrings(
            path, mcp_len, rev=rev, max_line_count=max_line_count
        )  # get MCP count dict
        mcp_cons_list = build_read_substring_cons_dict_list(
            mcp_count_dict, mcp_len
        )  # list of base conservation dicts for mcps
        cons_seq, cons_confs = build_cons_seq(
            mcp_cons_list,
            read_count,
            n_prop,
            do_not_include_list,
            max_line_count=max_line_count,
        )  # get list of max base conservations for each index
        # also get consensus sequence
        cons_seq_list.append(cons_seq)
        start_confs.append(np.mean(cons_confs))
        start_cons_lens.append(len(cons_seq))

    for i, end in enumerate(
        inf_point_list
    ):  # Looping through the post-inflection point mcps
        mcp_len = end + 5  # length of pre-inf mcps is inflection point + 5
        subs_len = start_cons_lens[i]  # length of respective pre-inf point sequence

        mcp_count_dict = fetch_read_substrings(
            path, subs_len, rev, mcp_len, max_line_count=max_line_count
        )
        mcp_cons_list = build_read_substring_cons_dict_list(mcp_count_dict, subs_len)
        cons_seq, cons_confs = build_cons_seq(
            mcp_cons_list,
            read_count,
            n_prop,
            do_not_include_list,
            subs_len,
            max_line_count=max_line_count,
        )

        end_confs.append(np.mean(cons_confs))

    diff_res = [
        start_confs[i] - end_confs[i] for i in range(len(start_confs))
    ]  # get differences between pre- and -post avg conservation values
    diff_res_sorted = sorted(
        diff_res, reverse=True
    )  # sort differences from highest to lowest

    ini_max_res = diff_res_sorted[0]  # maximum differences
    curr_max_index = diff_res.index(ini_max_res)  # index of maximum differences

    for res in diff_res_sorted[1:]:  # Loop through the rest of the differences
        curr_res_index = np.where(diff_res == res)[0][0]

        index_diff = inf_point_list[curr_max_index] - inf_point_list[curr_res_index]

        # if difference between the max and the current is negligible and the index of the current is earlier then..
        if ini_max_res - res < 0.05 and (index_diff <= 3 and index_diff > 0):
            curr_max_index = (
                curr_res_index  # replace the selected index with the current one
            )

    cutoff = (
        inf_point_list[curr_max_index] + 5
    )  # cutoff is the inflection point index + 5
    primer = cons_seq_list[
        curr_max_index
    ]  # grab the correct consensus sequence as primer

    # if the requested strand is reverse..
    if rev:
        primer = str(
            Seq(primer).complement()
        )  # ..get the complement of consensus sequence

    return cutoff, primer


def main():

    path, points, sample, output = parse_args()
    inf_df = pd.read_csv(points, sep="\t")

    f_slice = inf_df[inf_df.strand == "F"]  # get forward inflection points
    r_slice = inf_df[inf_df.strand == "R"]  # get reverse inflection points
    r_slice = r_slice.reset_index(drop=True)

    f_cutoff = ""
    r_cutoff = ""
    f_primer = ""
    r_primer = ""

    if not f_slice.empty:  # if there is a forward inflection point..
        inf_list = f_slice.inf_point.tolist()
        f_cutoff, f_primer = assess_inflection_point_mcp_for_sample(
            path, inf_list
        )  # .. assess and select

    if not r_slice.empty:  # if there is a reverse inflection point..
        inf_list = r_slice.inf_point.tolist()
        r_cutoff, r_primer = assess_inflection_point_mcp_for_sample(
            path, inf_list, rev=True
        )  # .. assess and select

    # Output cutoff point(s) to .txt file
    with open(f"{output}/{sample}_cutoff.txt", "w") as fw:
        if f_cutoff != "":
            fw.write(f"F: {f_cutoff}\n")
        if r_cutoff != "":
            fw.write(f"R: {r_cutoff}\n")

    # Output consensus primer sequence(s) to .fasta file
    with open(f"{output}/{sample}_auto_primers.fasta", "w") as fw:
        if f_cutoff != "":
            fw.write(f">F_auto\n{f_primer}\n")
        if r_cutoff != "":
            fw.write(f">R_auto\n{r_primer}\n")


if __name__ == "__main__":
    main()
