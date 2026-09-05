"""Input handling for pipette calibration CSV files."""

import csv
from typing import Dict, List, Optional

_REQUIRED_COLUMNS = [
    "PipetteID",
    "NominalVolume",
    "TargetVolume",
    "Weight_mg",
    "Temperature_C",
]


def read_calibration_data(filepath: str) -> List[Dict]:
    """Parse a CSV file and return a list of validated measurement dicts.

    Raises FileNotFoundError, ValueError, or KeyError on invalid data.
    """
    with open(filepath, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError("CSV file appears empty.")

        # Normalize headers: strip whitespace
        reader.fieldnames = [name.strip() for name in reader.fieldnames]

        # Check required columns
        missing = [col for col in _REQUIRED_COLUMNS if col not in reader.fieldnames]
        if missing:
            raise KeyError(f"Missing required columns: {missing}")

        rows = []
        for line_num, row in enumerate(reader, start=2):  # 1-indexed, header line 1
            try:
                # Strip whitespace from values
                cleaned = {k.strip(): v.strip() for k, v in row.items() if k is not None}
                # Convert and validate
                pip_id = str(cleaned["PipetteID"])
                if not pip_id:
                    raise ValueError("PipetteID cannot be empty")

                nominal = float(cleaned["NominalVolume"])
                target = float(cleaned["TargetVolume"])
                weight = float(cleaned["Weight_mg"])
                temp = float(cleaned["Temperature_C"])

                # Optional pressure
                pressure_str = cleaned.get("AirPressure_hPa", "").strip()
                if pressure_str:
                    pressure = float(pressure_str)
                else:
                    pressure = 1013.25

                if nominal <= 0 or target <= 0 or target > nominal:
                    raise ValueError(
                        f"Invalid nominal ({nominal}) or target ({target}) volume."
                    )
                if weight <= 0:
                    raise ValueError(f"Weight must be positive, got {weight}")
                if not (15.0 <= temp <= 30.0):
                    raise ValueError(f"Temperature {temp}°C out of valid range.")

                rows.append(
                    {
                        "PipetteID": pip_id,
                        "NominalVolume": nominal,
                        "TargetVolume": target,
                        "Weight_mg": weight,
                        "Temperature_C": temp,
                        "AirPressure_hPa": pressure,
                    }
                )
            except (ValueError, KeyError) as e:
                raise ValueError(f"Line {line_num}: {e}") from e

    if not rows:
        raise ValueError("No valid data rows found.")
    return rows
