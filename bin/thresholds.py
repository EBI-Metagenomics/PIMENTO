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

# used by bin.amplicon_utils.fetch_read_substrings
MAX_READ_COUNT = 300_000
# used by bin.amplicon_utils.build_cons_seq
CONSENSUS_BASE_THRESHOLD = 0.80

# used by pimento.standard_primer_strategy
MIN_STD_PRIMER_THRESHOLD = 0.60

# used by bin.standard_primer_matching.run_primer_matching_once
STD_PRIMER_READ_PREFIX_LENGTH = 50

# used by bin.are_there_primers.atp_in_this_sample
ATP_WINDOW_SIZE = 10
ATP_PREFIX_LENGTH = 100

# used by bin.generate_bcv.generate_bcv_for_single_strand
BCV_WINDOW_SIZE = 5

# used by bin.find_cutoffs.find_bcv_inflection_points
MIN_CUTOFF_WINDOW_START = 10
MAX_CUTOFF_WINDOW_START = 20
