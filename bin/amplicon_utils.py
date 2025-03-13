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

from collections import defaultdict, Counter
from pathlib import Path
import logging
import pyfastx

from bin.regex_ambiguous_bases import (
    AMBIGUOUS_BASES_DICT,
    AMBIGUOUS_BASES_DICT_REV,
)

logging.basicConfig(level=logging.DEBUG)


def get_read_count(read_path: Path, file_type: str = "fastq") -> int:
    """
    Get the read count of a FASTQ or FASTA file.

    :param read_path: The path to the FASTQ or FASTA file.
    :type read_path: Path
    :param fasta_type: The type of the file, either "fastq" or "fasta". Defaults to "fastq".
    :type fasta_type: str
    :return: The number of reads in the file.
    :rtype: int
    :raises ValueError: If the file type is not supported or the read count is not a positive integer.
    """
    read_count = 0

    if file_type == "fasta":
        fasta = pyfastx.Fasta(str(read_path), build_index=False)
        read_count = sum(1 for _ in fasta)
    elif file_type == "fastq":
        fastq = pyfastx.Fastq(str(read_path), build_index=False)
        read_count = sum(1 for _ in fastq)
    else:
        raise ValueError(
            f"Invalid file_type {file_type}, it needs to be either 'fasta' or 'fastq'"
        )

    if read_count <= 0:
        raise ValueError(f"Read count is not a positive integer: {read_count}")

    return read_count


def build_cons_seq(
    cons_list,
    read_count,
    cons_threshold=0.80,
    do_not_include=None,
    counter=1,
    max_line_count=0,
):
    """
    Generate consensus sequence using a list of base conservation dictionaries most likely
    generated by the `build_mcp_cons_dict_list()` function.
    Also returns a list containing the conservation value of the most conserved base at every
    position in the list of base conservation dictionaries.
    """

    cons_seq = ""
    cons_confs = []

    if do_not_include is None:
        do_not_include = []

    for count_dict in cons_list:
        max_count = 0
        cons_dict = defaultdict(float)

        if counter in do_not_include:
            counter += 1
            cons_seq += "N"
            continue

        for base, count in count_dict.items():
            if base not in ("A", "T", "C", "G"):
                continue

            if max_line_count == 0:
                cons_dict[base] = count / read_count
            else:
                cons_dict[base] = count / max_line_count

            if count > max_count:
                max_count = count

        counter += 1

        try:
            if max_line_count == 0:
                max_prop = max_count / read_count
            else:
                max_prop = max_count / max_line_count

            cons_bases = []
            curr_prop = 0.0
            sorted_cons_dict = dict(
                sorted(cons_dict.items(), key=lambda x: x[1], reverse=True)
            )

            for base, prop in sorted_cons_dict.items():
                cons_bases.append(base)
                curr_prop += prop
                if curr_prop >= cons_threshold:
                    break

            cons_bases = sorted(cons_bases)

            if len(cons_bases) == 1:
                cons_seq += cons_bases[0]
            else:
                amb_string = ",".join(cons_bases)
                amb_base = AMBIGUOUS_BASES_DICT_REV[amb_string]
                cons_seq += amb_base

        except ZeroDivisionError:
            max_prop = 0.0

        cons_confs.append(max_prop)

    return cons_seq, cons_confs


def primer_regex_query_builder(primer):
    """
    Takes an input nucleotide sequence that can contain IUPAC ambiguous codes
    Returns a string formatted as a regex query that considers the different
    potential bases valid at a position with am abiguity code.
    """

    query = ""

    for char in primer:
        if char in ("A", "C", "T", "G"):
            query += char
        else:
            query += str(AMBIGUOUS_BASES_DICT[char])

    query = f"(.*{query}){{e<=1}}"

    return query


def build_mcp_cons_dict_list(mcp_count_dict, mcp_len):
    """
    Generate list of dictionaries of base conservation for mcp output (mcp_cons_list)
    e.g. [{'A':0.9, 'C':0.1}, {'T':1.0}, ....] for every base position
    """

    mcp_cons_list = []

    for i in range(mcp_len):
        index_base_dict = defaultdict(int)
        for mcp in mcp_count_dict.keys():
            if len(mcp) < mcp_len:
                continue
            base = mcp[i]
            index_base_dict[base] += mcp_count_dict[mcp]
        mcp_cons_list.append(index_base_dict)

    return mcp_cons_list


def fetch_read_substrings(
    input_fastq: Path,
    prefix_len: int,
    rev: bool = False,
    start: int = 1,
    max_line_count: int = 0,
):
    """
    Generates the most common prefix sequences along with their counts in a fastq file.
    Outputs dictionary containing counts for each generated MCP in the fastq.
    """

    selected_lines = []

    fastq = pyfastx.Fastq(str(input_fastq), build_index=False)

    for read in fastq:
        sequence = read[1]  # the read sequence is the second element
        if not rev:
            selected_lines.append(sequence[start - 1 : start + prefix_len - 1])
        else:
            rev_sequence = sequence[::-1]
            selected_lines.append(rev_sequence[start - 1 : start + prefix_len - 1])
        if max_line_count != 0:
            if len(selected_lines) > max_line_count:
                break

    sequence_counts = Counter(selected_lines)
    mcp_count_dict = dict(
        sorted(sequence_counts.items(), key=lambda x: x[1], reverse=True)
    )

    return mcp_count_dict
