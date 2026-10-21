"""Rebuild the composite_solar cooling curve from Dere/Colgan/SPEX/cloudy_solar.

This is the canonical generator for the 151-point `t_composite` / `l_composite`
DATA block in src/physics/mod_radiative_cooling.t. Running this script with
`python3 tools/build_composite_solar.py` re-emits that DATA block on stdout.

Recipe
------
1. Source tables (log10 T, log10 Lambda per n_H^2, solar metallicity):
     - Dere_corona  (CHIANTI) : 101 pts, log T in [4.00, 9.00] step 0.05
     - Colgan       (LANL)    :  55 pts, log T in [4.065, 9.065]
     - SPEX         (Schure09): 110 pts, log T in [3.80, 8.16] step 0.04
     - cloudy_solar (Cloudy)  : 151 pts, log T in [1.0, 9.0] (= output grid)
2. Interpolate each source onto a fine grid `log_T_fine = arange(1.0, 9.001, 0.001)`.
3. Coronal weighted average (Dere:Colgan:SPEX = 6:4:2). Each weight ramps
   from 0 to its plateau value over the first 0.08 dex of its table's domain
   (smooth tanh ramp) and drops sharply to 0 outside the domain.
4. Zone blend (SPEX-only <-> coronal average):
       w_C = 0.5*(1 + tanh((log_T - 4.00) / 0.08))
5. Low-T blend (cloudy <-> zone blend):
       w_S = 0.5*(1 + tanh((log_T - 3.90) / 0.08))
6. Lower envelope: composite = max(composite, cloudy_solar).
7. Above log T = 4.48 the curve is forced to pure Dere_corona. A linear
   4.38 -> 4.48 ramp removes the join kink.
8. Savitzky-Golay smoothing, window_length = 201 fine points (0.2 dex),
   polyorder = 3.
9. Lower envelope with cloudy_solar reapplied after smoothing.
10. Sample the smoothed composite at the 151 cloudy_solar log-T nodes.

Output values are written to 3 decimal places (matching the historical
precision of the other cooling curves in the module).
"""
import numpy as np
from scipy.signal import savgol_filter


# --- Source tables (log10 T [K], log10 Lambda [erg cm^3 / s]) ---

t_Der = np.arange(4.00, 9.01, 0.05)
l_Der = np.array([-23.007, -22.554, -22.156, -21.833, -21.646, -21.616,
                  -21.684, -21.791, -21.876, -21.910, -21.900, -21.860,
                  -21.796, -21.717, -21.623, -21.521, -21.418, -21.331,
                  -21.277, -21.252, -21.242, -21.258, -21.279, -21.272,
                  -21.244, -21.215, -21.196, -21.183, -21.187, -21.241,
                  -21.352, -21.451, -21.487, -21.476, -21.452, -21.433,
                  -21.426, -21.420, -21.410, -21.401, -21.405, -21.423,
                  -21.450, -21.478, -21.511, -21.559, -21.635, -21.738,
                  -21.851, -21.955, -22.035, -22.090, -22.118, -22.124,
                  -22.116, -22.101, -22.083, -22.067, -22.057, -22.056,
                  -22.068, -22.101, -22.159, -22.240, -22.329, -22.408,
                  -22.471, -22.519, -22.550, -22.570, -22.580, -22.583,
                  -22.582, -22.575, -22.565, -22.554, -22.541, -22.527,
                  -22.514, -22.499, -22.485, -22.471, -22.456, -22.440,
                  -22.424, -22.407, -22.389, -22.371, -22.352, -22.332,
                  -22.312, -22.291, -22.269, -22.248, -22.226, -22.204,
                  -22.181, -22.159, -22.136, -22.113, -22.089])

t_Col = np.array([4.065, 4.142, 4.220, 4.298, 4.375, 4.453, 4.531, 4.608,
                  4.686, 4.764, 4.797, 4.830, 4.864, 4.897, 4.931, 4.964,
                  4.998, 5.031, 5.065, 5.176, 5.287, 5.398, 5.509, 5.620,
                  5.731, 5.842, 5.954, 6.065, 6.176, 6.287, 6.398, 6.509,
                  6.620, 6.731, 6.842, 6.954, 7.065, 7.176, 7.287, 7.398,
                  7.509, 7.620, 7.731, 7.842, 7.954, 8.065, 8.176, 8.287,
                  8.398, 8.509, 8.620, 8.731, 8.842, 8.954, 9.065])
l_Col = np.array([-22.189, -21.786, -21.604, -21.685, -21.764, -21.679,
                  -21.542, -21.380, -21.252, -21.176, -21.158, -21.145,
                  -21.135, -21.128, -21.125, -21.124, -21.126, -21.128,
                  -21.125, -21.090, -21.088, -21.195, -21.346, -21.348,
                  -21.317, -21.291, -21.289, -21.341, -21.431, -21.624,
                  -21.867, -22.029, -22.081, -22.061, -22.020, -22.000,
                  -22.052, -22.222, -22.415, -22.526, -22.569, -22.575,
                  -22.562, -22.540, -22.515, -22.489, -22.461, -22.429,
                  -22.394, -22.355, -22.313, -22.268, -22.222, -22.174,
                  -22.125])

t_Spx = np.arange(3.80, 8.17, 0.04)
l_Spx = np.array([-25.73, -25.04, -24.41, -23.83, -23.30, -22.82, -22.39,
                  -22.01, -21.68, -21.45, -21.32, -21.35, -21.43, -21.53,
                  -21.61, -21.66, -21.66, -21.59, -21.51, -21.41, -21.31,
                  -21.20, -21.11, -21.02, -20.94, -20.87, -20.82, -20.78,
                  -20.75, -20.75, -20.76, -20.78, -20.80, -20.80, -20.78,
                  -20.77, -20.76, -20.75, -20.75, -20.75, -20.79, -20.88,
                  -21.05, -21.23, -21.37, -21.46, -21.49, -21.51, -21.53,
                  -21.59, -21.65, -21.71, -21.74, -21.76, -21.77, -21.79,
                  -21.81, -21.83, -21.84, -21.85, -21.87, -21.91, -21.97,
                  -22.06, -22.15, -22.24, -22.31, -22.36, -22.40, -22.43,
                  -22.44, -22.44, -22.44, -22.43, -22.42, -22.42, -22.41,
                  -22.42, -22.43, -22.44, -22.46, -22.48, -22.51, -22.54,
                  -22.58, -22.62, -22.65, -22.67, -22.69, -22.70, -22.70,
                  -22.70, -22.70, -22.70, -22.69, -22.68, -22.67, -22.65,
                  -22.64, -22.63, -22.61, -22.60, -22.58, -22.57, -22.55,
                  -22.54, -22.52, -22.51, -22.49, -22.48])

t_cl = np.concatenate([np.arange(1.0, 8.001, 0.05),
                       np.arange(8.1, 9.01, 0.1)])
l_cl = np.array([-28.375, -28.251, -28.137, -28.029, -27.929, -27.834,
                 -27.745, -27.662, -27.584, -27.512, -27.445, -27.383,
                 -27.326, -27.273, -27.223, -27.175, -27.128, -27.079,
                 -27.027, -26.972, -26.911, -26.846, -26.777, -26.705,
                 -26.632, -26.554, -26.479, -26.407, -26.338, -26.274,
                 -26.213, -26.156, -26.101, -26.049, -25.999, -25.949,
                 -25.901, -25.852, -25.803, -25.754, -25.707, -25.662,
                 -25.621, -25.588, -25.561, -25.538, -25.518, -25.497,
                 -25.475, -25.452, -25.426, -25.400, -25.374, -25.333,
                 -25.295, -25.261, -25.228, -25.189, -25.136, -25.053,
                 -24.888, -24.454, -23.480, -22.562, -22.009, -21.826,
                 -21.840, -21.905, -21.956, -21.971, -21.958, -21.928,
                 -21.879, -21.810, -21.724, -21.623, -21.512, -21.404,
                 -21.321, -21.273, -21.250, -21.253, -21.275, -21.287,
                 -21.282, -21.275, -21.272, -21.267, -21.281, -21.357,
                 -21.496, -21.616, -21.677, -21.698, -21.708, -21.730,
                 -21.767, -21.793, -21.794, -21.787, -21.787, -21.802,
                 -21.826, -21.859, -21.911, -21.987, -22.082, -22.173,
                 -22.253, -22.325, -22.392, -22.448, -22.487, -22.512,
                 -22.524, -22.528, -22.524, -22.516, -22.507, -22.501,
                 -22.502, -22.511, -22.533, -22.565, -22.600, -22.630,
                 -22.648, -22.656, -22.658, -22.654, -22.647, -22.634,
                 -22.619, -22.602, -22.585, -22.566, -22.546, -22.525,
                 -22.505, -22.480, -22.465, -22.415, -22.365, -22.315,
                 -22.265, -22.215, -22.165, -22.115, -22.065, -22.015,
                 -21.965])


# --- Construction ---

log_T_fine = np.arange(1.0, 9.001, 0.001)
i_Der_f = np.interp(log_T_fine, t_Der, l_Der, left=np.nan, right=np.nan)
i_Col_f = np.interp(log_T_fine, t_Col, l_Col, left=np.nan, right=np.nan)
i_Spx_f = np.interp(log_T_fine, t_Spx, l_Spx, left=np.nan, right=np.nan)
i_cl_f  = np.interp(log_T_fine, t_cl,  l_cl)


def taper(log_T, lo, hi, w, ramp=0.08):
    """Per-curve weight profile: 0 outside [lo, hi], smooth tanh ramp over
    the first `ramp` dex above lo, plateau at w, sharp cutoff at hi."""
    out = np.zeros_like(log_T)
    m = (log_T >= lo) & (log_T <= hi)
    out[m] = w * 0.5 * (1 + np.tanh(8 * (log_T[m] - lo - ramp / 2) / ramp))
    out[log_T > lo + ramp] = np.where(log_T[log_T > lo + ramp] <= hi, w, 0)
    return out


w_D = taper(log_T_fine, 4.00,  9.00,  6.0)
w_C = taper(log_T_fine, 4.065, 9.065, 4.0)
w_S = taper(log_T_fine, 3.80,  8.16,  2.0)
w_D[np.isnan(i_Der_f)] = 0
w_C[np.isnan(i_Col_f)] = 0
w_S[np.isnan(i_Spx_f)] = 0

wt = w_D + w_C + w_S
coronal_avg = np.where(
    wt > 0,
    (w_D * np.nan_to_num(i_Der_f) +
     w_C * np.nan_to_num(i_Col_f) +
     w_S * np.nan_to_num(i_Spx_f)) / wt,
    np.nan,
)

spex_only = np.interp(log_T_fine, t_Spx, l_Spx, left=np.nan, right=np.nan)
w_C_blend = 0.5 * (1 + np.tanh((log_T_fine - 4.00) / 0.08))
zone_SC = np.where(
    ~np.isnan(coronal_avg),
    (1 - w_C_blend) * np.nan_to_num(spex_only, nan=-26.0) + w_C_blend * coronal_avg,
    np.nan_to_num(spex_only, nan=-26.0),
)

w_S_blend = 0.5 * (1 + np.tanh((log_T_fine - 3.90) / 0.08))
raw_blend = (1 - w_S_blend) * i_cl_f + w_S_blend * zone_SC

composite = np.maximum(raw_blend, i_cl_f)

hi = log_T_fine >= 4.48
composite[hi] = i_Der_f[hi]
for i in np.where((log_T_fine >= 4.38) & (log_T_fine < 4.48))[0]:
    f = (log_T_fine[i] - 4.38) / 0.10
    composite[i] = (1 - f) * composite[i] + f * i_Der_f[i]

composite_smooth = savgol_filter(composite, window_length=201, polyorder=3)
composite_smooth = np.maximum(composite_smooth, i_cl_f)

t_out = t_cl
l_out = np.interp(t_out, log_T_fine, composite_smooth)


# --- Emit Fortran DATA statements to stdout ---

def _emit(name, vals):
    print(f"  data    {name} /    &")
    for i in range(0, len(vals), 4):
        chunk = vals[i:i + 4]
        line = ", ".join(f"{v:9.3f}" for v in chunk)
        sep = ", &" if i + 4 < len(vals) else " / "
        print(f"     {line}   {sep}")
    print()


if __name__ == "__main__":
    print("  data n_composite / 151 /")
    print()
    _emit("t_composite", t_out)
    _emit("l_composite", l_out)
