"""Small HMI helpers required by the portable CEA preprocessing workflow."""

import numpy as np


def hmi_disambiguate_azimuth(azimuth, disambig, method: int = 2):
    """Apply the selected HMI disambiguation bit to an azimuth array."""

    method = 2 if method < 0 or method > 2 else int(method)
    azimuth = np.asarray(azimuth, dtype=float)
    disambig = np.nan_to_num(disambig, nan=0, posinf=0, neginf=0).astype(np.int32)
    if azimuth.shape != disambig.shape:
        raise ValueError("azimuth and disambig must have the same shape")

    corrected = azimuth.copy()
    corrected[((disambig >> method) & 1).astype(bool)] += 180.0
    return corrected
