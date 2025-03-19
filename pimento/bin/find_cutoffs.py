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

import numpy as np

from pimento.bin.thresholds import MAX_CUTOFF_WINDOW_START, MIN_CUTOFF_WINDOW_START


def find_bcv_inflection_points(bcv_df):
    """
    Find inflection points from an bcv_df file output by "generate_bcv.py"

    Takes the list of average base conservations and gets the derivative of the curve
    Keep any points of the derivative where value is above the 80th percentile

    Outputs a dictionary with key-val pairs where vals are lists:
        'strand' -> strand list
        'inf_point' -> inf_point list

    """

    inf_point_dict = defaultdict(list)
    start_indices = [int(i) for i in bcv_df.columns.tolist()]

    for i in range(len(bcv_df)):  # Loop through both possible strands of the bcv_df
        strand = bcv_df.index[i]
        props = bcv_df.iloc[i].tolist()
        props = [-val for val in props]

        prop_diff = np.diff(props) / np.diff(start_indices)  # Get the derivative
        infl_points = np.where(prop_diff > np.percentile(prop_diff, 80))[
            0
        ]  # Grab points above 80th percentile

        for ind in infl_points:
            inf_point = start_indices[ind]

            if (
                inf_point < MIN_CUTOFF_WINDOW_START
                or inf_point > MAX_CUTOFF_WINDOW_START
            ):  # Rule to facilitate results - won't accept
                continue  # points below index 10 or above index 20
                # 10 means a cutoff of 15 and 20 a cutoff of 25
                # literature points to no primers existing that are
                # shorter or bigger  than these lengths

            inf_point_dict["strand"].append(strand)
            inf_point_dict["inf_point"].append(inf_point)

    return inf_point_dict
