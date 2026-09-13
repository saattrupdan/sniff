"""Smoke coverage for the installed package and resource-backed CLI."""

import inspect
import json
import sys
import unittest
from contextlib import redirect_stdout
from importlib import resources
from io import StringIO
from unittest import mock

from sniff import analyze, catalogue, ptrms, viz


class _SingleArgumentTraversable:
    """Minimal Traversable fixture matching Python 3.9's joinpath API."""

    def __init__(self, *, children=None, content=None):
        self._children = children or {}
        self._content = content

    def joinpath(self, child):
        return self._children[child]

    def open(self, mode="r", encoding=None):
        if self._content is None:
            raise FileNotFoundError
        return StringIO(self._content)


class PackageSmokeTest(unittest.TestCase):
    def test_reference_resources_are_packaged(self):
        rate_constants = (
            resources.files("sniff")
            .joinpath("reference")
            .joinpath("rate_constants.json")
        )
        library = (
            resources.files("sniff")
            .joinpath("reference")
            .joinpath("ptrlibrary.csv")
        )
        compounds = (
            resources.files("sniff")
            .joinpath("reference")
            .joinpath("compound_catalogue.sqlite3")
        )
        self.assertTrue(rate_constants.is_file())
        self.assertTrue(library.is_file())
        self.assertTrue(compounds.is_file())
        self.assertGreater(
            int(catalogue.CompoundCatalogue().metadata()["species_count"]), 10_900
        )
        table = ptrms.load_rate_constants()
        self.assertIsNotNone(table)
        self.assertGreater(len(table["compounds"]), 100)

    def test_rate_constants_support_single_argument_traversable(self):
        table = {"compounds": [{"formula": "C2H6"}]}
        resource_tree = _SingleArgumentTraversable(
            children={
                "reference": _SingleArgumentTraversable(
                    children={
                        "rate_constants.json": _SingleArgumentTraversable(
                            content=json.dumps(table)
                        )
                    }
                )
            }
        )

        with mock.patch.object(ptrms.resources, "files", return_value=resource_tree):
            loaded = ptrms.load_rate_constants()

        self.assertEqual(loaded, table)

    def test_rates_command_needs_no_hdf5_fixture(self):
        output = StringIO()
        with (
            mock.patch.object(sys, "argv", ["sniff", "rates", "benzaldehyde"]),
            redirect_stdout(output),
        ):
            analyze.main()

        payload = json.loads(output.getvalue())
        self.assertIn("compounds", payload)
        self.assertTrue(payload["compounds"])

    def test_viz_waits_indefinitely_by_default(self):
        timeout_default = inspect.signature(viz.serve).parameters["timeout"].default
        self.assertIsNone(timeout_default)

        output = StringIO()
        with (
            mock.patch.object(sys, "argv", ["sniff", "viz", "--help"]),
            redirect_stdout(output),
            self.assertRaises(SystemExit) as raised,
        ):
            analyze.main()

        self.assertEqual(raised.exception.code, 0)
        self.assertIn("default: indefinitely", output.getvalue())

    def test_help_command_needs_no_hdf5_fixture(self):
        output = StringIO()
        with (
            mock.patch.object(sys, "argv", ["sniff", "--help"]),
            redirect_stdout(output),
            self.assertRaises(SystemExit) as raised,
        ):
            analyze.main()

        self.assertEqual(raised.exception.code, 0)
        self.assertIn("inspect", output.getvalue())


if __name__ == "__main__":
    unittest.main()
