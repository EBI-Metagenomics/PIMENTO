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
Unit tests for customizable threshold parameters in PIMENTO CLI commands.

Tests verify that user-supplied CLI parameters correctly override default
threshold values from pimento.bin.thresholds module.
"""

from pathlib import Path
from unittest.mock import patch, MagicMock

from click.testing import CliRunner

from pimento.pimento_cli import cli
from pimento.bin.thresholds import (
    MIN_STD_PRIMER_THRESHOLD,
    STD_PRIMER_READ_PREFIX_LENGTH,
    MAX_READ_COUNT,
    STD_PRIMER_ERROR_RATE,
)


class TestStandardPrimerThresholds:
    """Test threshold parameter customization for 'pimento std' command."""

    @patch("pimento.pimento_cli.write_std_output")
    @patch("pimento.pimento_cli.get_primer_props")
    @patch("pimento.pimento_cli.parse_std_primers")
    def test_std_minimum_primer_threshold_custom(
        self, mock_parse: MagicMock, mock_get_props: MagicMock, mock_write: MagicMock
    ) -> None:
        """
        Test that --minimum_primer_threshold parameter overrides MIN_STD_PRIMER_THRESHOLD default.

        Verifies that custom value (0.75) is passed to get_primer_props() instead of
        default value (0.60).

        :param mock_parse: Mock for parse_std_primers function.
        :type mock_parse: MagicMock
        :param mock_get_props: Mock for get_primer_props function.
        :type mock_get_props: MagicMock
        :param mock_write: Mock for write_std_output function.
        :type mock_write: MagicMock
        """
        # Setup mocks
        mock_parse.return_value = ({}, {}, 0)
        mock_get_props.return_value = []
        mock_write.return_value = (Path("test.fasta"), Path("test.txt"))

        runner = CliRunner()
        custom_threshold = 0.75

        result = runner.invoke(
            cli,
            [
                "std",
                "--input_fastq",
                "tests/fixtures/test.fastq.gz",
                "--output_prefix",
                "test_output",
                "--minimum_primer_threshold",
                str(custom_threshold),
                "--merged",
            ],
        )

        # Verify command executed successfully
        assert result.exit_code == 0

        # Verify get_primer_props was called with custom threshold
        assert mock_get_props.called
        call_args = mock_get_props.call_args[0]
        # Third positional argument is minimum_primer_threshold
        assert call_args[2] == custom_threshold
        assert call_args[2] != MIN_STD_PRIMER_THRESHOLD

    @patch("pimento.pimento_cli.write_std_output")
    @patch("pimento.pimento_cli.get_primer_props")
    @patch("pimento.pimento_cli.parse_std_primers")
    def test_std_primer_read_prefix_length_custom(
        self, mock_parse: MagicMock, mock_get_props: MagicMock, mock_write: MagicMock
    ) -> None:
        """
        Test that --std_primer_read_prefix_length parameter overrides STD_PRIMER_READ_PREFIX_LENGTH default.

        Verifies that custom value (75) is passed to get_primer_props() instead of
        default value (50).

        :param mock_parse: Mock for parse_std_primers function.
        :type mock_parse: MagicMock
        :param mock_get_props: Mock for get_primer_props function.
        :type mock_get_props: MagicMock
        :param mock_write: Mock for write_std_output function.
        :type mock_write: MagicMock
        """
        # Setup mocks
        mock_parse.return_value = ({}, {}, 0)
        mock_get_props.return_value = []
        mock_write.return_value = (Path("test.fasta"), Path("test.txt"))

        runner = CliRunner()
        custom_prefix_length = 75

        result = runner.invoke(
            cli,
            [
                "std",
                "--input_fastq",
                "tests/fixtures/test.fastq.gz",
                "--output_prefix",
                "test_output",
                "--std_primer_read_prefix_length",
                str(custom_prefix_length),
                "--merged",
            ],
        )

        # Verify command executed successfully
        assert result.exit_code == 0

        # Verify get_primer_props was called with custom prefix length
        assert mock_get_props.called
        call_args = mock_get_props.call_args[0]
        # Fourth positional argument is std_primer_read_prefix_length
        assert call_args[3] == custom_prefix_length
        assert call_args[3] != STD_PRIMER_READ_PREFIX_LENGTH

    @patch("pimento.pimento_cli.write_std_output")
    @patch("pimento.pimento_cli.get_primer_props")
    @patch("pimento.pimento_cli.parse_std_primers")
    def test_std_max_read_count_custom(
        self, mock_parse: MagicMock, mock_get_props: MagicMock, mock_write: MagicMock
    ) -> None:
        """
        Test that --max_read_count parameter overrides MAX_READ_COUNT default for std command.

        Verifies that custom value (100000) is passed to get_primer_props() instead of
        default value (300000).

        :param mock_parse: Mock for parse_std_primers function.
        :type mock_parse: MagicMock
        :param mock_get_props: Mock for get_primer_props function.
        :type mock_get_props: MagicMock
        :param mock_write: Mock for write_std_output function.
        :type mock_write: MagicMock
        """
        # Setup mocks
        mock_parse.return_value = ({}, {}, 0)
        mock_get_props.return_value = []
        mock_write.return_value = (Path("test.fasta"), Path("test.txt"))

        runner = CliRunner()
        custom_max_read_count = 100000

        result = runner.invoke(
            cli,
            [
                "std",
                "--input_fastq",
                "tests/fixtures/test.fastq.gz",
                "--output_prefix",
                "test_output",
                "--max_read_count",
                str(custom_max_read_count),
                "--merged",
            ],
        )

        # Verify command executed successfully
        assert result.exit_code == 0

        # Verify get_primer_props was called with custom max_read_count
        assert mock_get_props.called
        call_args = mock_get_props.call_args[0]
        # Fifth positional argument is max_read_count
        assert call_args[4] == custom_max_read_count
        assert call_args[4] != MAX_READ_COUNT

    @patch("pimento.pimento_cli.write_std_output")
    @patch("pimento.pimento_cli.get_primer_props")
    @patch("pimento.pimento_cli.parse_std_primers")
    def test_std_primer_error_rate_custom(
        self, mock_parse: MagicMock, mock_get_props: MagicMock, mock_write: MagicMock
    ) -> None:
        """
        Test that --std_primer_error_rate parameter overrides STD_PRIMER_ERROR_RATE default.

        Verifies that custom value (0.2) is passed to parse_std_primers() instead of
        default value (0.1).

        :param mock_parse: Mock for parse_std_primers function.
        :type mock_parse: MagicMock
        :param mock_get_props: Mock for get_primer_props function.
        :type mock_get_props: MagicMock
        :param mock_write: Mock for write_std_output function.
        :type mock_write: MagicMock
        """
        # Setup mocks
        mock_parse.return_value = ({}, {}, 0)
        mock_get_props.return_value = []
        mock_write.return_value = (Path("test.fasta"), Path("test.txt"))

        runner = CliRunner()
        custom_error_rate = 0.2

        result = runner.invoke(
            cli,
            [
                "std",
                "--input_fastq",
                "tests/fixtures/test.fastq.gz",
                "--output_prefix",
                "test_output",
                "--std_primer_error_rate",
                str(custom_error_rate),
                "--merged",
            ],
        )

        # Verify command executed successfully
        assert result.exit_code == 0

        # Verify parse_std_primers was called with custom error rate
        assert mock_parse.called
        call_args = mock_parse.call_args[0]
        # Second positional argument is std_primer_error_rate
        assert call_args[1] == custom_error_rate
        assert call_args[1] != STD_PRIMER_ERROR_RATE

    @patch("pimento.pimento_cli.write_std_output")
    @patch("pimento.pimento_cli.get_primer_props")
    @patch("pimento.pimento_cli.parse_std_primers")
    def test_std_all_custom_parameters(
        self, mock_parse: MagicMock, mock_get_props: MagicMock, mock_write: MagicMock
    ) -> None:
        """
        Test that all three threshold parameters work together with custom values.

        Verifies that --minimum_primer_threshold, --std_primer_read_prefix_length,
        and --max_read_count can all be customized simultaneously.

        :param mock_parse: Mock for parse_std_primers function.
        :type mock_parse: MagicMock
        :param mock_get_props: Mock for get_primer_props function.
        :type mock_get_props: MagicMock
        :param mock_write: Mock for write_std_output function.
        :type mock_write: MagicMock
        """
        # Setup mocks
        mock_parse.return_value = ({}, {}, 0)
        mock_get_props.return_value = []
        mock_write.return_value = (Path("test.fasta"), Path("test.txt"))

        runner = CliRunner()
        custom_threshold = 0.70
        custom_prefix_length = 80
        custom_max_read_count = 150000

        result = runner.invoke(
            cli,
            [
                "std",
                "--input_fastq",
                "tests/fixtures/test.fastq.gz",
                "--output_prefix",
                "test_output",
                "--minimum_primer_threshold",
                str(custom_threshold),
                "--std_primer_read_prefix_length",
                str(custom_prefix_length),
                "--max_read_count",
                str(custom_max_read_count),
                "--merged",
            ],
        )

        # Verify command executed successfully
        assert result.exit_code == 0

        # Verify get_primer_props was called with all custom values
        assert mock_get_props.called
        call_args = mock_get_props.call_args[0]

        assert call_args[2] == custom_threshold
        assert call_args[3] == custom_prefix_length
        assert call_args[4] == custom_max_read_count

        # Verify they differ from defaults
        assert call_args[2] != MIN_STD_PRIMER_THRESHOLD
        assert call_args[3] != STD_PRIMER_READ_PREFIX_LENGTH
        assert call_args[4] != MAX_READ_COUNT


class TestGenerateBCVThresholds:
    """Test threshold parameter customization for 'pimento gen_bcv' command."""

    @patch("click.Path.convert")
    @patch("pimento.pimento_cli.write_bcv_output")
    @patch("pimento.pimento_cli.generate_bcv_for_single_strand")
    def test_gen_bcv_max_read_count_custom(
        self,
        mock_generate_bcv: MagicMock,
        mock_write: MagicMock,
        mock_path_convert: MagicMock,
    ) -> None:
        """
        Test that --max_read_count parameter overrides MAX_READ_COUNT default for gen_bcv command.

        Verifies that custom value (200000) is passed to generate_bcv_for_single_strand()
        instead of default value (300000). Tests with strand='FR' which passes max_read_count.

        :param mock_generate_bcv: Mock for generate_bcv_for_single_strand function.
        :type mock_generate_bcv: MagicMock
        :param mock_write: Mock for write_bcv_output function.
        :type mock_write: MagicMock
        :param mock_path_convert: Mock for click.Path.convert to bypass file validation.
        :type mock_path_convert: MagicMock
        """
        # Mock path conversion to bypass file existence checks
        mock_path_convert.return_value = Path("mock_input.fastq.gz")

        # Setup mocks
        mock_generate_bcv.return_value = {}
        mock_write.return_value = MagicMock(to_csv=MagicMock())

        runner = CliRunner()
        custom_max_read_count = 200000

        result = runner.invoke(
            cli,
            [
                "gen_bcv",
                "--input_fastq",
                "mock_input.fastq.gz",
                "--strand",
                "FR",
                "--max_read_count",
                str(custom_max_read_count),
                "--output_prefix",
                "test_output",
            ],
        )

        # Verify command executed successfully
        assert result.exit_code == 0

        # Verify generate_bcv_for_single_strand was called with custom max_read_count
        # When strand='FR', the function is called twice (forward and reverse)
        assert mock_generate_bcv.called
        assert mock_generate_bcv.call_count == 2

        # Check both calls received custom max_read_count as second argument
        first_call_args = mock_generate_bcv.call_args_list[0][0]
        second_call_args = mock_generate_bcv.call_args_list[1][0]

        assert first_call_args[1] == custom_max_read_count
        assert second_call_args[1] == custom_max_read_count
        assert custom_max_read_count != MAX_READ_COUNT


class TestChoosePrimerCutoffThresholds:
    """Test threshold parameter customization for 'pimento choose_primer_cutoff' command."""

    @patch("click.Path.convert")
    @patch("pimento.pimento_cli.choose_cutoff_for_single_strand")
    @patch("pandas.read_csv")
    def test_choose_primer_cutoff_max_read_count_custom(
        self,
        mock_read_csv: MagicMock,
        mock_choose: MagicMock,
        mock_path_convert: MagicMock,
    ) -> None:
        """
        Test that --max_read_count parameter overrides MAX_READ_COUNT default for choose_primer_cutoff.

        Verifies that custom value (250000) is passed to choose_cutoff_for_single_strand()
        instead of default value (300000).

        :param mock_read_csv: Mock for pandas.read_csv function.
        :type mock_read_csv: MagicMock
        :param mock_choose: Mock for choose_cutoff_for_single_strand function.
        :type mock_choose: MagicMock
        :param mock_path_convert: Mock for click.Path.convert to bypass file validation.
        :type mock_path_convert: MagicMock
        """
        # Setup mocks - simulate dataframe with one forward primer cutoff
        import pandas as pd

        # Mock path conversion to bypass file existence checks
        mock_path_convert.return_value = Path("mock_file")

        mock_df = pd.DataFrame({"strand": ["F"], "inf_point": [15]})
        mock_read_csv.return_value = mock_df
        mock_choose.return_value = (15, "ATCGATCG")

        runner = CliRunner()
        custom_max_read_count = 250000

        result = runner.invoke(
            cli,
            [
                "choose_primer_cutoff",
                "--input_fastq",
                "mock_input.fastq.gz",
                "--primer_cutoffs",
                "mock_cutoffs.tsv",
                "--max_read_count",
                str(custom_max_read_count),
                "--output_prefix",
                "test_output",
            ],
        )

        # Verify command executed successfully
        assert result.exit_code == 0

        # Verify choose_cutoff_for_single_strand was called with custom max_read_count
        assert mock_choose.called
        call_args = mock_choose.call_args[0]
        # Third positional argument is max_read_count
        assert call_args[2] == custom_max_read_count
        assert call_args[2] != MAX_READ_COUNT


class TestAutoPipelineThresholds:
    """Test threshold parameter customization for 'pimento auto' command."""

    @patch("click.Path.convert")
    @patch("pimento.pimento_cli.pyfastx.Fasta")
    @patch("pimento.pimento_cli.choose_cutoff_for_single_strand")
    @patch("pimento.pimento_cli.find_bcv_inflection_points")
    @patch("pimento.pimento_cli.generate_bcv_for_single_strand")
    @patch("pandas.read_csv")
    def test_auto_max_read_count_custom(
        self,
        mock_read_csv: MagicMock,
        mock_generate_bcv: MagicMock,
        mock_find_cutoffs: MagicMock,
        mock_choose: MagicMock,
        mock_fasta: MagicMock,
        mock_path_convert: MagicMock,
    ) -> None:
        """
        Test that --max_read_count parameter propagates through auto command to subcommands.

        Verifies that custom value (175000) is passed through the entire auto pipeline
        (gen_bcv -> find_cutoffs -> choose_primer_cutoff).

        :param mock_read_csv: Mock for pandas.read_csv function.
        :type mock_read_csv: MagicMock
        :param mock_generate_bcv: Mock for generate_bcv_for_single_strand function.
        :type mock_generate_bcv: MagicMock
        :param mock_find_cutoffs: Mock for find_bcv_inflection_points function.
        :type mock_find_cutoffs: MagicMock
        :param mock_choose: Mock for choose_cutoff_for_single_strand function.
        :type mock_choose: MagicMock
        :param mock_fasta: Mock for pyfastx.Fasta class.
        :type mock_fasta: MagicMock
        :param mock_path_convert: Mock for click.Path.convert to bypass file validation.
        :type mock_path_convert: MagicMock
        """
        import pandas as pd

        # Mock path conversion to bypass file existence checks
        mock_path_convert.return_value = Path("mock_file")

        # Setup mocks
        mock_generate_bcv.return_value = {}
        mock_find_cutoffs.return_value = {"strand": ["F"], "inf_point": [15]}
        mock_choose.return_value = (15, "ATCGATCG")

        # Mock DataFrame for find_cutoffs
        mock_df = pd.DataFrame({"strand": ["F"], "inf_point": [15]})
        mock_read_csv.return_value = mock_df

        # Mock Fasta for final output reading
        mock_fasta_obj = MagicMock()
        mock_fasta_obj.__iter__ = MagicMock(return_value=iter([("F_auto", "ATCG")]))
        mock_fasta.return_value = mock_fasta_obj

        runner = CliRunner()
        custom_max_read_count = 175000

        result = runner.invoke(
            cli,
            [
                "auto",
                "--input_fastq",
                "mock_input.fastq.gz",
                "--strand",
                "FR",
                "--max_read_count",
                str(custom_max_read_count),
                "--output_prefix",
                "test_output",
            ],
        )

        # Verify command executed successfully
        assert result.exit_code == 0

        # Verify generate_bcv_for_single_strand was called with custom max_read_count
        # With strand='FR', it's called twice (forward and reverse)
        assert mock_generate_bcv.called
        assert mock_generate_bcv.call_count == 2

        # Check both calls have custom max_read_count
        first_gen_bcv_call = mock_generate_bcv.call_args_list[0][0]
        second_gen_bcv_call = mock_generate_bcv.call_args_list[1][0]
        assert first_gen_bcv_call[1] == custom_max_read_count
        assert second_gen_bcv_call[1] == custom_max_read_count

        # Verify choose_cutoff_for_single_strand was called with custom max_read_count
        assert mock_choose.called
        choose_call_args = mock_choose.call_args[0]
        assert choose_call_args[2] == custom_max_read_count
