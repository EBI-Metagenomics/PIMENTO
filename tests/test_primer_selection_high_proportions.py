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

"""
Unit tests for primer selection logic with significantly different proportions.

Tests verify that primers with significantly higher proportions (difference > 0.03)
are always selected, and that greedy length logic only applies when proportions
are similar (difference <= 0.03).
"""

from collections import defaultdict
from pathlib import Path
from unittest.mock import MagicMock, patch


class TestPrimerSelectionHighProportions:
    """Test primer selection logic handles high proportion primers correctly."""

    @patch("pimento.bin.standard_primer_matching.get_read_count")
    @patch("pimento.bin.standard_primer_matching.fetch_read_substrings")
    @patch("pimento.bin.standard_primer_matching.run_primer_matching_once")
    def test_high_proportion_primer_selected_longest_flag(
        self,
        mock_run_primer: MagicMock,
        mock_fetch: MagicMock,
        mock_read_count: MagicMock,
    ) -> None:
        """
        Test that a primer with significantly higher proportion is selected.

        When one primer has a proportion significantly higher than others
        (difference > 0.03), it should always be selected regardless of length
        or the greedy_primer_length_flag value.

        Primers:
        - primer_AF: 0.60 proportion (600 matches), 18bp
        - primer_BF: 0.90 proportion (900 matches), 25bp (HIGHEST - should be selected)
        - primer_CF: 0.65 proportion (650 matches), 20bp

        Flag: 'longest' (but shouldn't matter since 0.90 is clearly highest)
        Expected: primer_BF selected (0.90 proportion)
        """
        from pimento.bin.standard_primer_matching import get_primer_props

        # Mock the dependencies
        mock_read_count.return_value = 1000
        mock_fetch.return_value = {}  # Not used in this test

        # Mock primer matching results to return different counts
        # The function will be called once for each primer
        mock_run_primer.side_effect = [
            600,  # primer_AF: 0.60 proportion
            900,  # primer_BF: 0.90 proportion (highest)
            650,  # primer_CF: 0.65 proportion
        ]

        # Create primer dictionaries
        std_primer_dict_regex = defaultdict(dict)
        std_primer_dict = defaultdict(dict)

        std_primer_dict_regex["COI-5P"]["primer_AF"] = "A" * 18
        std_primer_dict_regex["COI-5P"]["primer_BF"] = "B" * 25
        std_primer_dict_regex["COI-5P"]["primer_CF"] = "C" * 20

        std_primer_dict["COI-5P"]["primer_AF"] = "A" * 18
        std_primer_dict["COI-5P"]["primer_BF"] = "B" * 25
        std_primer_dict["COI-5P"]["primer_CF"] = "C" * 20

        # Call get_primer_props
        result = get_primer_props(
            std_primer_dict_regex=std_primer_dict_regex,
            std_primer_dict=std_primer_dict,
            input_fastq=Path("test.fastq.gz"),
            min_std_primer_threshold=0.50,
            std_primer_read_prefix_length=50,
            max_read_count=1000,
            greedy_primer_length_flag="longest",
            merged=True,
            threads=1,
        )

        # Verify primer_BF (highest proportion) was selected
        # Result format: [region_name, {primer_name: proportion}]
        assert len(result) == 2
        assert result[0] == "COI-5P"
        assert "primer_BF" in result[1]
        assert result[1]["primer_BF"] == 0.90

    @patch("pimento.bin.standard_primer_matching.get_read_count")
    @patch("pimento.bin.standard_primer_matching.fetch_read_substrings")
    @patch("pimento.bin.standard_primer_matching.run_primer_matching_once")
    def test_high_proportion_primer_selected_shortest_flag(
        self,
        mock_run_primer: MagicMock,
        mock_fetch: MagicMock,
        mock_read_count: MagicMock,
    ) -> None:
        """
        Test that high proportion overrides shortest flag.

        Even with greedy_primer_length_flag='shortest', a primer with
        significantly higher proportion should be selected.

        Primers:
        - primer_AF: 0.60 proportion (600 matches), 18bp (shortest)
        - primer_BF: 0.90 proportion (900 matches), 25bp (HIGHEST prop - should be selected)
        - primer_CF: 0.65 proportion (650 matches), 20bp

        Flag: 'shortest' (but 0.90 is so much higher it should override)
        Expected: primer_BF selected (0.90 proportion)
        """
        from pimento.bin.standard_primer_matching import get_primer_props

        # Mock the dependencies
        mock_read_count.return_value = 1000
        mock_fetch.return_value = {}

        # Mock primer matching results
        mock_run_primer.side_effect = [600, 900, 650]

        # Create primer dictionaries
        std_primer_dict_regex = defaultdict(dict)
        std_primer_dict = defaultdict(dict)

        std_primer_dict_regex["COI-5P"]["primer_AF"] = "A" * 18
        std_primer_dict_regex["COI-5P"]["primer_BF"] = "B" * 25
        std_primer_dict_regex["COI-5P"]["primer_CF"] = "C" * 20

        std_primer_dict["COI-5P"]["primer_AF"] = "A" * 18
        std_primer_dict["COI-5P"]["primer_BF"] = "B" * 25
        std_primer_dict["COI-5P"]["primer_CF"] = "C" * 20

        # Call get_primer_props with 'shortest' flag
        result = get_primer_props(
            std_primer_dict_regex=std_primer_dict_regex,
            std_primer_dict=std_primer_dict,
            input_fastq=Path("test.fastq.gz"),
            min_std_primer_threshold=0.50,
            std_primer_read_prefix_length=50,
            max_read_count=1000,
            greedy_primer_length_flag="shortest",
            merged=True,
            threads=1,
        )

        # Verify primer_BF (highest proportion) was selected despite 'shortest' flag
        assert len(result) == 2
        assert result[0] == "COI-5P"
        assert "primer_BF" in result[1]
        assert result[1]["primer_BF"] == 0.90

    @patch("pimento.bin.standard_primer_matching.get_read_count")
    @patch("pimento.bin.standard_primer_matching.fetch_read_substrings")
    @patch("pimento.bin.standard_primer_matching.run_primer_matching_once")
    def test_high_proportion_late_in_iteration(
        self,
        mock_run_primer: MagicMock,
        mock_fetch: MagicMock,
        mock_read_count: MagicMock,
    ) -> None:
        """
        Test that high proportion primer is selected even if it appears last.

        This specifically tests the bug where primers appearing later in iteration
        with higher proportions were incorrectly skipped.

        Primers (in iteration order):
        - primer_AF: 0.60 proportion (600 matches)
        - primer_BF: 0.62 proportion (620 matches)
        - primer_CF: 0.95 proportion (950 matches) (appears last, HIGHEST - should be selected)

        Expected: primer_CF selected despite appearing last
        """
        from pimento.bin.standard_primer_matching import get_primer_props

        # Mock the dependencies
        mock_read_count.return_value = 1000
        mock_fetch.return_value = {}

        # Mock primer matching results
        mock_run_primer.side_effect = [600, 620, 950]

        # Create primer dictionaries
        std_primer_dict_regex = defaultdict(dict)
        std_primer_dict = defaultdict(dict)

        std_primer_dict_regex["COI-5P"]["primer_AF"] = "A" * 20
        std_primer_dict_regex["COI-5P"]["primer_BF"] = "B" * 20
        std_primer_dict_regex["COI-5P"]["primer_CF"] = "C" * 20

        std_primer_dict["COI-5P"]["primer_AF"] = "A" * 20
        std_primer_dict["COI-5P"]["primer_BF"] = "B" * 20
        std_primer_dict["COI-5P"]["primer_CF"] = "C" * 20

        # Call get_primer_props
        result = get_primer_props(
            std_primer_dict_regex=std_primer_dict_regex,
            std_primer_dict=std_primer_dict,
            input_fastq=Path("test.fastq.gz"),
            min_std_primer_threshold=0.50,
            std_primer_read_prefix_length=50,
            max_read_count=1000,
            greedy_primer_length_flag="longest",
            merged=True,
            threads=1,
        )

        # Verify primer_CF (highest proportion) was selected
        assert len(result) == 2
        assert result[0] == "COI-5P"
        assert "primer_CF" in result[1]
        assert result[1]["primer_CF"] == 0.95

    @patch("pimento.bin.standard_primer_matching.get_read_count")
    @patch("pimento.bin.standard_primer_matching.fetch_read_substrings")
    @patch("pimento.bin.standard_primer_matching.run_primer_matching_once")
    def test_greedy_logic_still_works_for_close_proportions(
        self,
        mock_run_primer: MagicMock,
        mock_fetch: MagicMock,
        mock_read_count: MagicMock,
    ) -> None:
        """
        Test that greedy length logic still applies when proportions are close.

        When primers have similar proportions (difference <= 0.03), the greedy
        length flag should determine which primer is selected.

        Primers:
        - primer_AF: 0.60 proportion (600 matches), 18bp (shorter)
        - primer_BF: 0.62 proportion (620 matches), 25bp (longer, diff=0.02 <= 0.03)

        Flag: 'longest'
        Expected: primer_BF selected (0.62 proportion, longer)
        """
        from pimento.bin.standard_primer_matching import get_primer_props

        # Mock the dependencies
        mock_read_count.return_value = 1000
        mock_fetch.return_value = {}

        # Mock primer matching results - proportions differ by 0.02 (< 0.03)
        mock_run_primer.side_effect = [600, 620]

        # Create primer dictionaries with different lengths
        std_primer_dict_regex = defaultdict(dict)
        std_primer_dict = defaultdict(dict)

        std_primer_dict_regex["COI-5P"]["primer_AF"] = "A" * 18
        std_primer_dict_regex["COI-5P"]["primer_BF"] = "B" * 25

        std_primer_dict["COI-5P"]["primer_AF"] = "A" * 18
        std_primer_dict["COI-5P"]["primer_BF"] = "B" * 25

        # Call get_primer_props with 'longest' flag
        result = get_primer_props(
            std_primer_dict_regex=std_primer_dict_regex,
            std_primer_dict=std_primer_dict,
            input_fastq=Path("test.fastq.gz"),
            min_std_primer_threshold=0.50,
            std_primer_read_prefix_length=50,
            max_read_count=1000,
            greedy_primer_length_flag="longest",
            merged=True,
            threads=1,
        )

        # Verify primer_BF (longer) was selected
        assert len(result) == 2
        assert result[0] == "COI-5P"
        assert "primer_BF" in result[1]
        assert result[1]["primer_BF"] == 0.62

    @patch("pimento.bin.standard_primer_matching.get_read_count")
    @patch("pimento.bin.standard_primer_matching.fetch_read_substrings")
    @patch("pimento.bin.standard_primer_matching.run_primer_matching_once")
    def test_equal_proportions_greedy_applies(
        self,
        mock_run_primer: MagicMock,
        mock_fetch: MagicMock,
        mock_read_count: MagicMock,
    ) -> None:
        """
        Test that greedy logic applies when proportions are exactly equal.

        When primers have identical proportions, the greedy length flag
        should determine selection.

        Primers:
        - primer_AF: 0.60 proportion (600 matches), 18bp (shorter)
        - primer_BF: 0.60 proportion (600 matches), 25bp (longer)

        Flag: 'longest'
        Expected: primer_BF selected (equal proportion, but longer)
        """
        from pimento.bin.standard_primer_matching import get_primer_props

        # Mock the dependencies
        mock_read_count.return_value = 1000
        mock_fetch.return_value = {}

        # Mock primer matching results - identical proportions
        mock_run_primer.side_effect = [600, 600]

        # Create primer dictionaries with different lengths
        std_primer_dict_regex = defaultdict(dict)
        std_primer_dict = defaultdict(dict)

        std_primer_dict_regex["COI-5P"]["primer_AF"] = "A" * 18
        std_primer_dict_regex["COI-5P"]["primer_BF"] = "B" * 25

        std_primer_dict["COI-5P"]["primer_AF"] = "A" * 18
        std_primer_dict["COI-5P"]["primer_BF"] = "B" * 25

        # Call get_primer_props with 'longest' flag
        result = get_primer_props(
            std_primer_dict_regex=std_primer_dict_regex,
            std_primer_dict=std_primer_dict,
            input_fastq=Path("test.fastq.gz"),
            min_std_primer_threshold=0.50,
            std_primer_read_prefix_length=50,
            max_read_count=1000,
            greedy_primer_length_flag="longest",
            merged=True,
            threads=1,
        )

        # Verify primer_BF (longer) was selected
        assert len(result) == 2
        assert result[0] == "COI-5P"
        assert "primer_BF" in result[1]
        assert result[1]["primer_BF"] == 0.60
