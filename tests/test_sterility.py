"""Unit tests for sterility test calculator module."""

import pytest

from lab_tools.sterility.sterility import (
    calculate_sample_size,
    calculate_volume_per_container,
    calculate_direct_inoculation_volume,
    main,
)


class TestSampleSize:
    """Tests for calculate_sample_size."""

    # ---- Liquid / USP ----
    def test_liquid_usp_small_batch(self):
        size, warn = calculate_sample_size(50, "liquid", "USP")
        assert size == 5
        assert warn is None

    def test_liquid_usp_tiny_batch(self):
        size, warn = calculate_sample_size(2, "liquid", "USP")
        assert size == 2
        assert warn is not None
        assert "exceeds batch size" in warn

    def test_liquid_usp_mid_range(self):
        size, _ = calculate_sample_size(250, "liquid", "USP")
        assert size == 10

    def test_liquid_usp_upper_boundary(self):
        size, _ = calculate_sample_size(500, "liquid", "USP")
        assert size == 10

    def test_liquid_usp_large(self):
        size, _ = calculate_sample_size(600, "liquid", "USP")
        assert size == 12  # 2% of 600 = 12

    def test_liquid_usp_very_large_capped(self):
        size, _ = calculate_sample_size(2000, "liquid", "USP")
        assert size == 20  # 2% of 2000 = 40 → capped at 20

    # ---- Solid / USP ----
    def test_solid_usp_small(self):
        size, _ = calculate_sample_size(50, "solid", "USP")
        assert size == 2  # 5% of 50 = 2.5 → max(2,2)=2

    def test_solid_usp_boundary(self):
        size, _ = calculate_sample_size(200, "solid", "USP")
        assert size == 10  # 5% of 200 = 10, max(2,10)=10

    def test_solid_usp_large(self):
        size, _ = calculate_sample_size(300, "solid", "USP")
        assert size == 10

    # ---- Ophthalmic / USP (same as solid) ----
    def test_ophthalmic_usp_small(self):
        size, _ = calculate_sample_size(60, "ophthalmic", "USP")
        assert size == 3  # 5% of 60 = 3, max(2,3)=3

    def test_ophthalmic_usp_large(self):
        size, _ = calculate_sample_size(250, "ophthalmic", "USP")
        assert size == 10

    # ---- EP (identical rules in this implementation) ----
    def test_ep_liquid(self):
        size, _ = calculate_sample_size(80, "liquid", "EP")
        assert size == 8  # 10% of 80 = 8, max(4,8)=8

    # ---- Edge Cases ----
    def test_batch_size_zero(self):
        with pytest.raises(ValueError, match="positive integer"):
            calculate_sample_size(0, "liquid")

    def test_batch_size_negative(self):
        with pytest.raises(ValueError):
            calculate_sample_size(-5, "liquid")

    def test_invalid_product_type(self):
        with pytest.raises(ValueError, match="Invalid product type"):
            calculate_sample_size(100, "gel")

    def test_invalid_pharmacopeia(self):
        with pytest.raises(ValueError, match="Pharmacopeia"):
            calculate_sample_size(100, "liquid", "EPX")

    def test_batch_size_not_int(self):
        with pytest.raises(ValueError):
            calculate_sample_size(10.5, "liquid")  # type: ignore


class TestVolumePerContainer:
    """Tests for calculate_volume_per_container."""

    def test_liquid_small_container(self):
        vol, desc = calculate_volume_per_container("liquid", 50)
        assert vol == 50
        assert "entire contents" in desc

    def test_liquid_large_container(self):
        vol, desc = calculate_volume_per_container("liquid", 250)
        assert vol == 10.0
        assert "10 mL" in desc

    def test_ophthalmic_small_container(self):
        vol, desc = calculate_volume_per_container("ophthalmic", 4)
        assert vol == 4
        assert "entire" in desc

    def test_ophthalmic_large_container(self):
        vol, _ = calculate_volume_per_container("ophthalmic", 10)
        assert vol == 5.0

    def test_solid_not_allowed(self):
        with pytest.raises(ValueError, match="only applicable"):
            calculate_volume_per_container("solid", 100)

    def test_negative_volume(self):
        with pytest.raises(ValueError):
            calculate_volume_per_container("liquid", -5)


class TestDirectInoculation:
    """Tests for calculate_direct_inoculation_volume and direct mode CLI."""

    def test_direct_below_one_ml(self):
        vol, desc = calculate_direct_inoculation_volume(0.5)
        assert vol == 0.5
        assert "entire" in desc

    def test_direct_one_to_five_half(self):
        vol, desc = calculate_direct_inoculation_volume(4)
        assert vol == 2.0
        assert "half" in desc

    def test_direct_five_to_twenty(self):
        vol, desc = calculate_direct_inoculation_volume(10)
        assert vol == 2.0
        assert "2 mL" in desc

    def test_direct_twenty_or_more(self):
        vol, desc = calculate_direct_inoculation_volume(100)
        assert vol == 10.0
        assert "10 mL" in desc

    def test_direct_invalid_volume(self):
        with pytest.raises(ValueError):
            calculate_direct_inoculation_volume(0)

    def test_cli_direct_reports_volume(self, capsys):
        rc = main([
            "size", "--method", "direct", "--fill-volume", "100",
        ])
        assert rc == 0
        out = capsys.readouterr().out
        assert "direct inoculation" in out
        assert "10.0 mL" in out

    def test_cli_direct_with_batch_reports_units(self, capsys):
        rc = main([
            "--method", "direct", "--fill-volume", "100",
            "--batch-size", "100", "--product-type", "liquid",
        ])
        assert rc == 0
        out = capsys.readouterr().out
        assert "Minimum number of units to test: 10" in out
        assert "10.0 mL" in out
