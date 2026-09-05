import pytest
from lab_tools.microbiology.mpn import compute_mpn, format_result


class TestTypicalCases:
    def test_5tube_mixed(self):
        """Standard 5-tube 10/1/0.1 mL with 5-2-0 pattern approx 49 MPN."""
        res = compute_mpn(5, [10.0, 1.0, 0.1], [5, 2, 0])
        assert res["method"] == "mle"
        assert res["mpn"] is not None
        # Expect ~49, loose tolerance
        assert 40.0 <= res["mpn"] <= 60.0
        assert res["lower_ci"] > 0.0
        assert res["mpn"] <= res["upper_ci"]
        assert res["warning"] is None

    def test_3tube_mixed(self):
        """3-tube 10/1/0.1 mL with 3-1-0."""
        res = compute_mpn(3, [10.0, 1.0, 0.1], [3, 1, 0])
        assert res["method"] == "mle"
        # Typical table value ~43, accept wide range
        assert 20.0 <= res["mpn"] <= 80.0
        assert res["lower_ci"] < res["mpn"] < res["upper_ci"]

    def test_unit_factor(self):
        """Check that per_unit works (per 1 mL)"""
        res100 = compute_mpn(5, [10.0, 1.0, 0.1], [5, 2, 0], per_unit=100.0)
        res1   = compute_mpn(5, [10.0, 1.0, 0.1], [5, 2, 0], per_unit=1.0)
        assert abs(res100["mpn"] - res1["mpn"] * 100.0) < 0.01


class TestBoundaryCases:
    def test_all_negative(self):
        """5 tubes all zero -> MPN=0, upper CI computed."""
        res = compute_mpn(5, [10.0, 1.0, 0.1], [0, 0, 0])
        assert res["mpn"] == 0.0
        assert res["lower_ci"] == 0.0
        # Upper CI approx 5.4 per 100 mL
        assert 4.5 <= res["upper_ci"] <= 6.0
        assert res["warning"] is not None

    def test_all_positive(self):
        """All tubes positive -> MPN > measurable, lower CI given."""
        res = compute_mpn(5, [10.0, 1.0, 0.1], [5, 5, 5])
        assert res["mpn"] is None
        assert res["upper_ci"] is None
        assert res["lower_ci"] is not None and res["lower_ci"] > 0.0
        assert "All tubes positive" in res["warning"]


class TestInputValidation:
    def test_length_mismatch(self):
        with pytest.raises(ValueError, match="same length"):
            compute_mpn(5, [10.0, 1.0], [5, 2, 0])

    def test_positive_out_of_range(self):
        with pytest.raises(ValueError, match="out of range"):
            compute_mpn(5, [10.0, 1.0, 0.1], [6, 0, 0])

    def test_negative_positive_count(self):
        with pytest.raises(ValueError):
            compute_mpn(5, [10.0], [-1])

    def test_zero_or_negative_volume(self):
        with pytest.raises(ValueError, match="positive, got"):
            compute_mpn(5, [10.0, 0.0, 0.1], [5, 0, 0])

    def test_invalid_tubes(self):
        with pytest.raises(ValueError):
            compute_mpn(0, [10.0], [0])


class TestFormatting:
    def test_format_normal(self):
        res = compute_mpn(5, [10.0, 1.0, 0.1], [5, 2, 0])
        out = format_result(res)
        assert "MPN" in out
        assert "95% CI" in out
        assert "MLE" in out

    def test_format_warning(self):
        res = compute_mpn(5, [10.0, 1.0, 0.1], [0, 0, 0])
        out = format_result(res)
        assert "Warning:" in out
