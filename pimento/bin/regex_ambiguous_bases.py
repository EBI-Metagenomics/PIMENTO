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

AMBIGUOUS_BASES_DICT = {
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACTG]",
}

AMBIGUOUS_BASES_DICT_REV = {
    "A,G": "R",
    "C,T": "Y",
    "C,G": "S",
    "A,T": "W",
    "G,T": "K",
    "A,C": "M",
    "C,G,T": "B",
    "A,G,T": "D",
    "A,C,T": "H",
    "A,C,G": "V",
    "A,C,G,T": "N",
}
