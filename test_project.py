import unittest

import numpy as np
import pandas as pd

from analysis import carrier_fraction
from common import RESULTS


class CarrierCoverageTests(unittest.TestCase):
    def test_closed_form_and_bounds(self):
        self.assertAlmostEqual(carrier_fraction(0.4), 0.64)
        self.assertAlmostEqual(carrier_fraction(0.4, f=0.1), 0.616)
        self.assertEqual(carrier_fraction(0), 0)
        self.assertEqual(carrier_fraction(1), 1)
        with self.assertRaises(ValueError):
            carrier_fraction(1.1)

    def test_published_scenario_table(self):
        scenarios = pd.read_csv(RESULTS / "scenarios.csv")
        self.assertEqual(len(scenarios), 36)
        evaluated = scenarios.apply(lambda row: carrier_fraction(row.allele_mass, row.F), axis=1)
        np.testing.assert_allclose(scenarios.carrier_coverage, evaluated, atol=1e-6)
