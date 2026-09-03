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
Unit tests for primer input type handling in PIMENTO CLI.

Tests verify that the primer input parameter correctly handles both:
1. Directory input (default behavior)
2. FASTA file input (new feature)
3. Invalid file extensions (error handling)
"""

from pathlib import Path
from unittest.mock import patch, MagicMock

from click.testing import CliRunner

from pimento.pimento_cli import cli


class TestPrimerInputType:
    """Test primer input type handling for 'pimento std' command."""

    @patch("pimento.pimento_cli.write_std_output")
    @patch("pimento.pimento_cli.get_primer_props")
    @patch("pimento.pimento_cli.parse_std_primers")
    def test_primer_input_as_fasta_file(
        self, mock_parse: MagicMock, mock_get_props: MagicMock, mock_write: MagicMock
    ) -> None:
        """
        Test that a FASTA file can be used as primer input via -p parameter.

        Verifies that when a FASTA file is provided as the primer_input,
        the parse_std_primers function receives a file path (not a directory).

        :param mock_parse: Mock for parse_std_primers function.
        :type mock_parse: MagicMock
        :param mock_get_props: Mock for get_primer_props function.
        :type mock_get_props: MagicMock
        :param mock_write: Mock for write_std_output function.
        :type mock_write: MagicMock
        """
        # Setup mocks
        mock_parse.return_value = ({}, {}, 2)
        mock_get_props.return_value = []
        mock_write.return_value = (Path("test.fasta"), Path("test.txt"))

        runner = CliRunner()
        primer_fasta_file = Path("pimento/standard_primers/V3-V4.fasta")

        result = runner.invoke(
            cli,
            [
                "std",
                "--input_fastq",
                "tests/fixtures/test.fastq.gz",
                "--primer_input",
                str(primer_fasta_file),
                "--output_prefix",
                "test_output",
                "--merged",
            ],
        )

        # Verify command executed successfully
        assert result.exit_code == 0

        # Verify parse_std_primers was called with a file path
        assert mock_parse.called
        call_args = mock_parse.call_args[0]
        primer_input_arg = call_args[0]

        # Verify that the primer_input is a file (not a directory)
        assert primer_input_arg.is_file(), f"Expected file, got: {primer_input_arg}"
        assert (
            primer_input_arg == primer_fasta_file
        ), f"Expected {primer_fasta_file}, got {primer_input_arg}"

    @patch("pimento.pimento_cli.write_std_output")
    @patch("pimento.pimento_cli.get_primer_props")
    @patch("pimento.pimento_cli.parse_std_primers")
    def test_primer_input_as_directory(
        self, mock_parse: MagicMock, mock_get_props: MagicMock, mock_write: MagicMock
    ) -> None:
        """
        Test that directory input works (default behavior).

        Verifies that when no -p parameter is provided, the default primer directory
        is used and parse_std_primers receives a directory path.

        :param mock_parse: Mock for parse_std_primers function.
        :type mock_parse: MagicMock
        :param mock_get_props: Mock for get_primer_props function.
        :type mock_get_props: MagicMock
        :param mock_write: Mock for write_std_output function.
        :type mock_write: MagicMock
        """
        # Setup mocks
        mock_parse.return_value = ({}, {}, 10)
        mock_get_props.return_value = []
        mock_write.return_value = (Path("test.fasta"), Path("test.txt"))

        runner = CliRunner()

        # Don't provide --primer_input, so it uses the default directory
        result = runner.invoke(
            cli,
            [
                "std",
                "--input_fastq",
                "tests/fixtures/test.fastq.gz",
                "--output_prefix",
                "test_output",
                "--merged",
            ],
        )

        # Verify command executed successfully
        assert result.exit_code == 0

        # Verify parse_std_primers was called
        assert mock_parse.called
        call_args = mock_parse.call_args[0]
        primer_input_arg = call_args[0]

        # Verify that the primer_input is a directory (default behavior)
        assert primer_input_arg.is_dir(), f"Expected directory, got: {primer_input_arg}"

    def test_primer_input_invalid_file_extension(self, tmp_path: Path) -> None:
        """
        Test that invalid file extensions are rejected with proper error message.

        Verifies that when a file with an unsupported extension (not .fa, .fna, .fasta,
        or their .gz variants) is provided, a ValueError is raised with the expected message.

        :param tmp_path: Pytest fixture providing temporary directory.
        :type tmp_path: Path
        """
        # Create a temporary file with invalid extension
        invalid_file = tmp_path / "primers.txt"
        invalid_file.write_text(">primer1\nATCGATCG\n")

        runner = CliRunner()

        result = runner.invoke(
            cli,
            [
                "std",
                "--input_fastq",
                "tests/fixtures/test.fastq.gz",
                "--primer_input",
                str(invalid_file),
                "--output_prefix",
                "test_output",
                "--merged",
            ],
        )

        # Verify command failed
        assert result.exit_code != 0

        # Verify error message contains expected text
        # The error is in the exception, not the output
        assert result.exception is not None
        assert "is neither a directory nor a fasta file" in str(result.exception)

    def test_primer_input_fasta_gz_extension(self, tmp_path: Path) -> None:
        """
        Test that .fasta.gz file extension is accepted.

        Verifies that compressed FASTA files (.fasta.gz) are correctly recognized
        as valid primer input files.

        :param tmp_path: Pytest fixture providing temporary directory.
        :type tmp_path: Path
        """
        import gzip

        # Create a temporary .fasta.gz file
        fasta_gz_file = tmp_path / "primers.fasta.gz"
        with gzip.open(fasta_gz_file, "wt") as f:
            f.write(">341F\nCCTACGGGNGGCWGCAG\n")

        # We need to mock since we're not actually using this file for matching
        with (
            patch("pimento.pimento_cli.write_std_output") as mock_write,
            patch("pimento.pimento_cli.get_primer_props") as mock_get_props,
            patch("pimento.pimento_cli.parse_std_primers") as mock_parse,
        ):
            mock_parse.return_value = ({}, {}, 1)
            mock_get_props.return_value = []
            mock_write.return_value = (Path("test.fasta"), Path("test.txt"))

            runner = CliRunner()

            result = runner.invoke(
                cli,
                [
                    "std",
                    "--input_fastq",
                    "tests/fixtures/test.fastq.gz",
                    "--primer_input",
                    str(fasta_gz_file),
                    "--output_prefix",
                    "test_output",
                    "--merged",
                ],
            )

            # Verify command executed successfully
            assert result.exit_code == 0

            # Verify parse_std_primers was called with the .fasta.gz file
            assert mock_parse.called
            call_args = mock_parse.call_args[0]
            primer_input_arg = call_args[0]
            assert primer_input_arg.is_file()
            assert primer_input_arg.suffix == ".gz"
            assert str(primer_input_arg).endswith(".fasta.gz")

    def test_primer_count_matches_fasta_sequences(self) -> None:
        """
        Test that primer count from parse_std_primers matches the number of sequences in input FASTA.

        This test directly calls parse_std_primers with the V3-V4.fasta file and verifies
        that the returned primer_count matches the actual number of sequences in the file.
        """
        from pimento.bin.standard_primer_matching import parse_std_primers

        primer_fasta_file = Path("pimento/standard_primers/V3-V4.fasta")

        # V3-V4.fasta contains 8 primer sequences (3 forward, 5 reverse)
        expected_primer_count = 8

        # Call parse_std_primers directly
        std_primer_dict_regex, std_primer_dict, primer_count = parse_std_primers(
            primer_fasta_file, std_primer_error_rate=0.1, merged=False
        )

        # Verify the primer_count matches expected
        assert (
            primer_count == expected_primer_count
        ), f"Expected {expected_primer_count} primers, but got {primer_count}"

        # Additional verification: check that dictionaries are not empty
        assert (
            len(std_primer_dict_regex) > 0
        ), "std_primer_dict_regex should not be empty"
        assert len(std_primer_dict) > 0, "std_primer_dict should not be empty"
