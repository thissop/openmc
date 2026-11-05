import math
import unittest

from metrics import (
    figures_of_merit,
    median_positive_finite,
    normalized_pattern,
    relative_errors,
    rms_fractional_bias,
)


class MetricsTests(unittest.TestCase):
    def test_zero_mean_is_unresolved(self):
        errors = relative_errors([2.0, 0.0], [0.5, 0.0])
        self.assertAlmostEqual(errors[0], 0.25)
        self.assertTrue(math.isinf(errors[1]))

    def test_figure_of_merit(self):
        self.assertEqual(figures_of_merit([0.5, math.inf], 4.0), [1.0, 0.0])

    def test_pattern_normalization(self):
        self.assertEqual(normalized_pattern([1.0, 3.0]), [0.25, 0.75])
        self.assertTrue(all(math.isnan(x) for x in normalized_pattern([0.0, 0.0])))

    def test_fractional_bias(self):
        self.assertAlmostEqual(rms_fractional_bias([2.0, 4.0], [1.0, 2.0]), 1.0)

    def test_positive_finite_median(self):
        self.assertEqual(
            median_positive_finite([None, 0.0, 1.0, 3.0, math.inf]), 2.0
        )


if __name__ == "__main__":
    unittest.main()
