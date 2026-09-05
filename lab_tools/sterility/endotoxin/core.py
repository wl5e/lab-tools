"""Core bacterial endotoxin limit and MVD calculations per USP <85>."""

from typing import Optional, Dict


def calc_endotoxin_limit_product(
    dose_mg_per_kg_hour: float,
    concentration_mg_per_ml: Optional[float] = None,
    k: float = 5.0,
) -> Dict[str, float]:
    """
    Calculate endotoxin limit for a drug product given maximum mass dose.

    Parameters
    ----------
    dose_mg_per_kg_hour : float
        Maximum dose of active substance in mg per kg body weight per hour.
    concentration_mg_per_ml : float, optional
        Concentration of active substance in mg/mL. If provided, returns limit in EU/mL;
        otherwise returns limit in EU/mg of active.
    k : float, default 5.0
        K factor in EU/kg (e.g., 5 for IV, 0.2 for intrathecal). Must be > 0.

    Returns
    -------
    dict with keys:
        'endotoxin_limit'  : limit in EU/mg (if no concentration) or EU/mL (if concentration given)
        'unit'             : 'EU/mg' or 'EU/mL'
        'K'                : used K value
        'dose_mg_per_kg_hour' : input dose
        'concentration_mg_per_ml' : concentration if provided, else None

    Raises
    ------
    ValueError for invalid inputs.
    """
    if dose_mg_per_kg_hour <= 0:
        raise ValueError("Maximum dose (mg/kg/h) must be positive.")
    if k <= 0:
        raise ValueError("K factor must be positive.")
    if concentration_mg_per_ml is not None and concentration_mg_per_ml <= 0:
        raise ValueError("Concentration must be positive if provided.")

    # EL per mg active = K / dose (mg/kg/h)
    el_per_mg = k / dose_mg_per_kg_hour

    if concentration_mg_per_ml is None:
        return {
            "endotoxin_limit": round(el_per_mg, 6),
            "unit": "EU/mg",
            "K": k,
            "dose_mg_per_kg_hour": dose_mg_per_kg_hour,
            "concentration_mg_per_ml": None,
        }
    else:
        # limit (EU/mL) = el_per_mg * concentration_mg_per_ml
        el_per_ml = el_per_mg * concentration_mg_per_ml
        return {
            "endotoxin_limit": round(el_per_ml, 6),
            "unit": "EU/mL",
            "K": k,
            "dose_mg_per_kg_hour": dose_mg_per_kg_hour,
            "concentration_mg_per_ml": concentration_mg_per_ml,
        }


def calc_endotoxin_limit_volume(
    dose_ml_per_kg_hour: float,
    k: float = 5.0,
) -> Dict[str, float]:
    """
    Calculate endotoxin limit directly from dose volume rate.

    Parameters
    ----------
    dose_ml_per_kg_hour : float
        Maximum volume of product in mL per kg body weight per hour.
    k : float, default 5.0
        K factor in EU/kg.

    Returns
    -------
    dict with endotoxin_limit (EU/mL), unit, etc.
    """
    if dose_ml_per_kg_hour <= 0:
        raise ValueError("Dose volume (mL/kg/h) must be positive.")
    if k <= 0:
        raise ValueError("K factor must be positive.")

    el_per_ml = k / dose_ml_per_kg_hour
    return {
        "endotoxin_limit": round(el_per_ml, 6),
        "unit": "EU/mL",
        "K": k,
        "dose_ml_per_kg_hour": dose_ml_per_kg_hour,
    }


def calc_mvd(
    endotoxin_limit_eu_per_ml: float,
    lysate_sensitivity_eu_per_ml: float,
) -> float:
    """
    Calculate Maximum Valid Dilution for the BET assay.

    MVD = endotoxin_limit / lysate_sensitivity

    Parameters
    ----------
    endotoxin_limit_eu_per_ml : float
        Endotoxin limit of the product in EU/mL.
    lysate_sensitivity_eu_per_ml : float
        Labeled sensitivity of the lysate reagent in EU/mL.

    Returns
    -------
    float: MVD (dimensionless dilution factor).

    Raises
    ------
    ValueError if inputs are non-positive or invalid.
    """
    if endotoxin_limit_eu_per_ml <= 0:
        raise ValueError("Endotoxin limit must be positive.")
    if lysate_sensitivity_eu_per_ml <= 0:
        raise ValueError("Lysate sensitivity must be positive.")

    mvd = endotoxin_limit_eu_per_ml / lysate_sensitivity_eu_per_ml
    return mvd
