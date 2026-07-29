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
Unit tests for --greedy_primer_length_flag parameter in PIMENTO CLI.

Tests verify that the greedy_primer_length_flag parameter correctly selects
either the longest or shortest primer when multiple primers are viable candidates
(i.e., when the difference in proportions is small, <= 0.03).
"""

from pathlib import Path

from click.testing import CliRunner

from pimento.pimento_cli import cli


class TestGreedyPrimerLengthFlag:
    """Test greedy_primer_length_flag parameter for 'pimento std' command."""

    def test_greedy_longest_choice1(self, tmp_path: Path) -> None:
        """
        Test that --greedy_primer_length_flag longest selects longer primers.

        Uses test_longer_primer_std_choice1.fastq.gz which contains reads matching
        both GHalF (23bp) and COL6 (26bp) forward primers with similar proportions.
        When using 'longest' flag, COL6 should be selected as the forward primer.

        Expected output:
        - Forward: COL6|COI-5P|F (26bp) - TYTCHACAAAYCATAAAGAYATYGG
        - Reverse: MLepF1-Rev|COI-5P|R (22bp) - CGTGGAAAWGCTATATCWGGTG
        """
        runner = CliRunner()
        output_prefix = tmp_path / "test_output"

        result = runner.invoke(
            cli,
            [
                "std",
                "--input_fastq",
                "tests/fixtures/test_longer_primer_std_choice1.fastq.gz",
                "--primers_dir",
                "tests/fixtures/greedy_std_primer_len_choice_primers/",
                "--output_prefix",
                str(output_prefix),
                "--merged",
                "--greedy_primer_length_flag",
                "longest",
            ],
        )

        # Verify command executed successfully
        assert result.exit_code == 0

        # Read output FASTA file
        fasta_path = Path(f"{output_prefix}_std_primers.fasta")
        assert fasta_path.exists()

        # Parse FASTA file
        fasta_records = {}
        with open(fasta_path) as f:
            current_header = None
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    current_header = line[1:]  # Remove '>'
                    fasta_records[current_header] = ""
                elif current_header:
                    fasta_records[current_header] += line

        # Verify exactly 2 primers are found
        assert len(fasta_records) == 2

        # Verify COL6 (longer forward primer) is selected
        assert "COL6|COI-5P|F" in fasta_records
        assert fasta_records["COL6|COI-5P|F"] == "TYTCHACAAAYCATAAAGAYATYGG"

        # Verify MLepF1-Rev (reverse primer) is selected
        assert "MLepF1-Rev|COI-5P|R" in fasta_records
        assert fasta_records["MLepF1-Rev|COI-5P|R"] == "CGTGGAAAWGCTATATCWGGTG"

        # Verify GHalF (shorter forward primer) is NOT selected
        assert "GHalF|COI-5P|F" not in fasta_records

    def test_greedy_shortest_choice1(self, tmp_path: Path) -> None:
        """
        Test that --greedy_primer_length_flag shortest selects shorter primers.

        Uses test_longer_primer_std_choice1.fastq.gz which contains reads matching
        both GHalF (23bp) and COL6 (26bp) forward primers with similar proportions.
        When using 'shortest' flag, GHalF should be selected as the forward primer.

        Expected output:
        - Forward: GHalF|COI-5P|F (23bp) - TCAACAAATCATAAAGATATYGG
        - Reverse: MLepF1-Rev|COI-5P|R (22bp) - CGTGGAAAWGCTATATCWGGTG
        """
        runner = CliRunner()
        output_prefix = tmp_path / "test_output"

        result = runner.invoke(
            cli,
            [
                "std",
                "--input_fastq",
                "tests/fixtures/test_longer_primer_std_choice1.fastq.gz",
                "--primers_dir",
                "tests/fixtures/greedy_std_primer_len_choice_primers/",
                "--output_prefix",
                str(output_prefix),
                "--merged",
                "--greedy_primer_length_flag",
                "shortest",
            ],
        )

        # Verify command executed successfully
        assert result.exit_code == 0

        # Read output FASTA file
        fasta_path = Path(f"{output_prefix}_std_primers.fasta")
        assert fasta_path.exists()

        # Parse FASTA file
        fasta_records = {}
        with open(fasta_path) as f:
            current_header = None
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    current_header = line[1:]  # Remove '>'
                    fasta_records[current_header] = ""
                elif current_header:
                    fasta_records[current_header] += line

        # Verify exactly 2 primers are found
        assert len(fasta_records) == 2

        # Verify GHalF (shorter forward primer) is selected
        assert "GHalF|COI-5P|F" in fasta_records
        assert fasta_records["GHalF|COI-5P|F"] == "TCAACAAATCATAAAGATATYGG"

        # Verify MLepF1-Rev (reverse primer) is selected
        assert "MLepF1-Rev|COI-5P|R" in fasta_records
        assert fasta_records["MLepF1-Rev|COI-5P|R"] == "CGTGGAAAWGCTATATCWGGTG"

        # Verify COL6 (longer forward primer) is NOT selected
        assert "COL6|COI-5P|F" not in fasta_records

    def test_greedy_longest_choice2(self, tmp_path: Path) -> None:
        """
        Test that --greedy_primer_length_flag longest selects longer primers for both strands.

        Uses test_longer_primer_std_choice2.fastq.gz which contains reads matching:
        - Forward: M13F (18bp) and BirdF1_t1 (43bp) with similar proportions
        - Reverse: Fish_R1-T1 (38bp) and VR1_t1 (44bp) with similar proportions
        When using 'longest' flag, both longer primers should be selected.

        Expected output:
        - Forward: BirdF1_t1|COI-5P|F (43bp) - TGTAAAACGACGGCCAGTTCTCCAACCACAAAGACATTGGCAC
        - Reverse: VR1_t1|COI-5P|R (44bp) - CAGGAAACAGCTATGACTAGACTTCTGGGTGGCCAAAGAATCA
        """
        runner = CliRunner()
        output_prefix = tmp_path / "test_output"

        result = runner.invoke(
            cli,
            [
                "std",
                "--input_fastq",
                "tests/fixtures/test_longer_primer_std_choice2.fastq.gz",
                "--primers_dir",
                "tests/fixtures/greedy_std_primer_len_choice_primers/",
                "--output_prefix",
                str(output_prefix),
                "--merged",
                "--greedy_primer_length_flag",
                "longest",
            ],
        )

        # Verify command executed successfully
        assert result.exit_code == 0

        # Read output FASTA file
        fasta_path = Path(f"{output_prefix}_std_primers.fasta")
        assert fasta_path.exists()

        # Parse FASTA file
        fasta_records = {}
        with open(fasta_path) as f:
            current_header = None
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    current_header = line[1:]  # Remove '>'
                    fasta_records[current_header] = ""
                elif current_header:
                    fasta_records[current_header] += line

        # Verify exactly 2 primers are found
        assert len(fasta_records) == 2

        # Verify BirdF1_t1 (longer forward primer) is selected
        assert "BirdF1_t1|COI-5P|F" in fasta_records
        assert (
            fasta_records["BirdF1_t1|COI-5P|F"]
            == "TGTAAAACGACGGCCAGTTCTCCAACCACAAAGACATTGGCAC"
        )

        # Verify VR1_t1 (longer reverse primer) is selected
        assert "VR1_t1|COI-5P|R" in fasta_records
        assert (
            fasta_records["VR1_t1|COI-5P|R"]
            == "CAGGAAACAGCTATGACTAGACTTCTGGGTGGCCAAAGAATCA"
        )

        # Verify M13F (shorter forward primer) is NOT selected
        assert "M13F|COI-5P|F" not in fasta_records

        # Verify Fish_R1-T1 (shorter reverse primer) is NOT selected
        assert "Fish_R1-T1|COI-5P|R" not in fasta_records

    def test_greedy_shortest_choice2(self, tmp_path: Path) -> None:
        """
        Test that --greedy_primer_length_flag shortest selects shorter primers for both strands.

        Uses test_longer_primer_std_choice2.fastq.gz which contains reads matching:
        - Forward: M13F (18bp) and BirdF1_t1 (43bp) with similar proportions
        - Reverse: Fish_R1-T1 (38bp) and VR1_t1 (44bp) with similar proportions
        When using 'shortest' flag, both shorter primers should be selected.

        Expected output:
        - Forward: M13F|COI-5P|F (18bp) - TGTAAAACGACGGCCAGT
        - Reverse: Fish_R1-T1|COI-5P|R (38bp) - CAGGAAACAGCTATGACTAGACTTCTGGGTGGCCAAAG
        """
        runner = CliRunner()
        output_prefix = tmp_path / "test_output"

        result = runner.invoke(
            cli,
            [
                "std",
                "--input_fastq",
                "tests/fixtures/test_longer_primer_std_choice2.fastq.gz",
                "--primers_dir",
                "tests/fixtures/greedy_std_primer_len_choice_primers/",
                "--output_prefix",
                str(output_prefix),
                "--merged",
                "--greedy_primer_length_flag",
                "shortest",
            ],
        )

        # Verify command executed successfully
        assert result.exit_code == 0

        # Read output FASTA file
        fasta_path = Path(f"{output_prefix}_std_primers.fasta")
        assert fasta_path.exists()

        # Parse FASTA file
        fasta_records = {}
        with open(fasta_path) as f:
            current_header = None
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    current_header = line[1:]  # Remove '>'
                    fasta_records[current_header] = ""
                elif current_header:
                    fasta_records[current_header] += line

        # Verify exactly 2 primers are found
        assert len(fasta_records) == 2

        # Verify M13F (shorter forward primer) is selected
        assert "M13F|COI-5P|F" in fasta_records
        assert fasta_records["M13F|COI-5P|F"] == "TGTAAAACGACGGCCAGT"

        # Verify Fish_R1-T1 (shorter reverse primer) is selected
        assert "Fish_R1-T1|COI-5P|R" in fasta_records
        assert (
            fasta_records["Fish_R1-T1|COI-5P|R"]
            == "CAGGAAACAGCTATGACTAGACTTCTGGGTGGCCAAAG"
        )

        # Verify BirdF1_t1 (longer forward primer) is NOT selected
        assert "BirdF1_t1|COI-5P|F" not in fasta_records

        # Verify VR1_t1 (longer reverse primer) is NOT selected
        assert "VR1_t1|COI-5P|R" not in fasta_records
