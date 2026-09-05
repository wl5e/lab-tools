"""Core calculations for gravimetric pipette calibration."""

import math
from typing import Dict, List, Tuple

# ISO 8655-2:2002 tolerances for adjustable pipettes (systematic / random error in %)
# Values are maximum permissible error for *as-found* testing.
_TOLERANCE_TABLE = [
    (0.1, 1, 8.0, 8.0),
    (1, 2, 3.0, 3.0),
    (2, 5, 2.0, 2.0),
    (5, 10, 1.5, 1.5),
    (10, 20, 1.0, 0.8),
    (20, 100, 0.8, 0.3),
    (100, 200, 0.6, 0.3),
    (200, 1000, 0.6, 0.2),
    (1000, 10000, 0.5, 0.13),
]


def _water_density(temp_c: float) -> float:
    """Return water density in g/mL at given temperature (°C).

    Formula from ISO 8655-6:2002 Annex B, valid for 15–30 °C.
    """
    if not (15.0 <= temp_c <= 30.0):
        raise ValueError(f"Temperature {temp_c}°C is outside the valid range 15–30°C.")
    rho = (
        0.999849
        + 6.856e-5 * temp_c
        - 8.441e-6 * temp_c**2
        + 5.241e-8 * temp_c**3
        - 8.467e-10 * temp_c**4
    )
    return rho


def get_z_factor(temp_c: float, pressure_hpa: float = 1013.25) -> float:
    """Calculate the gravimetric conversion factor Z.

    Z converts mass (mg) directly to volume (µL) for water.
    Implements air buoyancy correction using standard brass weight density.
    """
    rho_water = _water_density(temp_c)  # g/mL
    # air density (g/mL) at given temp and pressure, simplified
    rho_air = 0.0012 * (pressure_hpa / 1013.25) * (293.15 / (temp_c + 273.15))
    rho_brass = 8.0  # g/mL

    # Z = volume_in_mL_per_gram * 1000 to get µL/mg? Actually:
    # V_mL = mass_g * (1 - rho_air/rho_brass) / (rho_water - rho_air)
    # mass_mg = mass_g * 1000, so V_µL = mass_mg * Z where Z = (1 - rho_air/rho_brass) / (rho_water - rho_air) / 1000
    # Wait: V_mL = mass_g * Z_mL_g  where Z_mL_g = (1 - rho_air/rho_brass) / (rho_water - rho_air)
    # Then V_µL = V_mL * 1000 = mass_g * Z_mL_g * 1000 = mass_mg * Z_mL_g / 1000 * 1000? Let's derive:
    # mass_g = mass_mg / 1000, so V_µL = (mass_mg/1000) * Z_mL_g * 1000 = mass_mg * Z_mL_g.
    # Therefore Z = Z_mL_g.
    z = (1.0 - rho_air / rho_brass) / (rho_water - rho_air)
    return z


def mass_to_volume(mass_mg: float, temp_c: float, pressure_hpa: float = 1013.25) -> float:
    """Convert mass (mg) to volume (µL) using the Z factor."""
    if mass_mg <= 0:
        raise ValueError("Mass must be positive.")
    return mass_mg * get_z_factor(temp_c, pressure_hpa)


def get_tolerances(nominal_ul: float) -> Tuple[float, float]:
    """Return (max_accuracy_pct, max_cv_pct) for a given nominal volume (µL)."""
    for low, high, acc, cv in _TOLERANCE_TABLE:
        if low <= nominal_ul < high:
            return acc, cv
    raise ValueError(f"Nominal volume {nominal_ul} µL out of supported range (0.1–10000 µL).")


def _compute_group_stats(
    target: float, measurements: List[Dict]
) -> Dict:
    """Compute statistics for a single target volume group."""
    volumes = [
        mass_to_volume(m["Weight_mg"], m["Temperature_C"], m.get("AirPressure_hPa", 1013.25))
        for m in measurements
    ]
    n = len(volumes)
    if n < 3:
        raise ValueError(
            f"Need at least 3 measurements per target point, got {n} for target {target} µL."
        )

    mean_vol = sum(volumes) / n
    if mean_vol == 0:
        raise ValueError("Mean volume is zero; cannot compute CV.")
    # systematic error (accuracy) in %
    accuracy_pct = (mean_vol - target) / target * 100.0
    # random error (CV) in %
    variance = sum((v - mean_vol) ** 2 for v in volumes) / (n - 1)
    std_dev = math.sqrt(variance)
    cv_pct = (std_dev / mean_vol) * 100.0

    return {
        "MeanVolume_ul": round(mean_vol, 3),
        "Accuracy_pct": round(accuracy_pct, 2),
        "CV_pct": round(cv_pct, 2),
        "N": n,
    }


def analyze_calibration(raw_data: List[Dict]) -> List[Dict]:
    """Analyze a list of measurement dicts and return a results table.

    Expects each dict to have keys: PipetteID, NominalVolume, TargetVolume,
    Weight_mg, Temperature_C, and optionally AirPressure_hPa.
    """
    # Group by (pipette_id, nominal_volume, target_volume)
    groups: Dict[Tuple[str, float, float], List[Dict]] = {}
    for row in raw_data:
        key = (row["PipetteID"], float(row["NominalVolume"]), float(row["TargetVolume"]))
        groups.setdefault(key, []).append(row)

    results = []
    for (pip_id, nominal, target), measurements in sorted(groups.items()):
        try:
            stats = _compute_group_stats(target, measurements)
        except ValueError as e:
            # Mark as fail with error message
            results.append(
                {
                    "PipetteID": pip_id,
                    "NominalVolume": nominal,
                    "TargetVolume": target,
                    "N": len(measurements),
                    "MeanVolume_ul": "N/A",
                    "Accuracy_pct": "Error",
                    "CV_pct": "Error",
                    "PassFail": f"FAIL ({str(e)})",
                }
            )
            continue

        # Determine pass/fail based on tolerances for the nominal volume
        try:
            max_acc_pct, max_cv_pct = get_tolerances(nominal)
        except ValueError:
            max_acc_pct, max_cv_pct = (0.0, 0.0)  # will fail
        passed = abs(stats["Accuracy_pct"]) <= max_acc_pct and stats["CV_pct"] <= max_cv_pct
        status = "PASS" if passed else "FAIL"

        results.append(
            {
                "PipetteID": pip_id,
                "NominalVolume": nominal,
                "TargetVolume": target,
                "N": stats["N"],
                "MeanVolume_ul": stats["MeanVolume_ul"],
                "Accuracy_pct": stats["Accuracy_pct"],
                "CV_pct": stats["CV_pct"],
                "PassFail": status,
            }
        )

    # Sort for consistent output
    results.sort(key=lambda r: (r["PipetteID"], r["NominalVolume"], r["TargetVolume"]))
    return results
