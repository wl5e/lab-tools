"""Tests for endotoxin_calc core module."""

import pytest

from lab_tools.sterility.endotoxin.core import (
    calc_endotoxin_limit_product,
    calc_endotoxin_limit_volume,
    calc_mvd,
)


class TestEndotoxinLimitProduct:
    def test_basic_with_concentration(self):
        result = calc_endotoxin_limit_product(
            dose_mg_per_kg_hour=10.0, concentration_mg_per_ml=2.0, k=5.0
        )
        assert result["unit"] == "EU/mL"
        assert result["endotoxin_limit"] == pytest.approx(1.0)  # K=5, dose=10 => el_per_mg=0.5, *2=1.0
        assert result["K"] == 5.0

    def test_basic_without_concentration(self):
        result = calc_endotoxin_limit_product(dose_mg_per_kg_hour=10.0, k=5.0)
        assert result["unit"] == "EU/mg"
        assert result["endotoxin_limit"] == pytest.approx(0.5)
        assert result["concentration_mg_per_ml"] is None

    def test_intrathecal_low_k(self):
        result = calc_endotoxin_limit_product(
            dose_mg_per_kg_hour=5.0, concentration_mg_per_ml=1.0, k=0.2
        )
        # el_per_mg = 0.2/5 = 0.04, *1 = 0.04 EU/mL
        assert result["endotoxin_limit"] == pytest.approx(0.04)

    def test_zero_dose_raises(self):
        with pytest.raises(ValueError, match="Maximum dose.*positive"):
            calc_endotoxin_limit_product(dose_mg_per_kg_hour=0.0)

    def test_negative_dose_raises(self):
        with pytest.raises(ValueError, match="Maximum dose.*positive"):
            calc_endotoxin_limit_product(dose_mg_per_kg_hour=-5.0)

    def test_negative_k_raises(self):
        with pytest.raises(ValueError, match="K factor.*positive"):
            calc_endotoxin_limit_product(dose_mg_per_kg_hour=10.0, k=-1.0)

    def test_zero_concentration_raises(self):
        with pytest.raises(ValueError, match="Concentration.*positive"):
            calc_endotoxin_limit_product(
                dose_mg_per_kg_hour=10.0, concentration_mg_per_ml=0.0
            )


class TestEndotoxinLimitVolume:
    def test_typical(self):
        result = calc_endotoxin_limit_volume(dose_ml_per_kg_hour=5.0, k=5.0)
        assert result["unit"] == "EU/mL"
        assert result["endotoxin_limit"] == pytest.approx(1.0)
        assert result["K"] == 5.0

    def test_small_volume(self):
        result = calc_endotoxin_limit_volume(dose_ml_per_kg_hour=0.5, k=0.2)
        assert result["endotoxin_limit"] == pytest.approx(0.4)  # 0.2/0.5

    def test_zero_volume_raises(self):
        with pytest.raises(ValueError, match="Dose volume.*positive"):
            calc_endotoxin_limit_volume(0.0)

    def test_negative_k_raises(self):
        with pytest.raises(ValueError, match="K factor.*positive"):
            calc_endotoxin_limit_volume(4.0, k=0.0)


class TestMVD:
    def test_normal(self):
        assert calc_mvd(0.5, 0.125) == pytest.approx(4.0)

    def test_small_limit(self):
        assert calc_mvd(0.005, 0.03) == pytest.approx(0.166666, abs=1e-4)

    def test_zero_limit_raises(self):
        with pytest.raises(ValueError, match="Endotoxin limit.*positive"):
            calc_mvd(0.0, 0.125)

    def test_negative_sensitivity_raises(self):
        with pytest.raises(ValueError, match="Lysate sensitivity.*positive"):
            calc_mvd(0.5, -0.125)
