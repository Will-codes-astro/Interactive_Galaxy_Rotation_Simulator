import urllib.request
import urllib.error
import ssl
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp

# ==============================================================================
# ARCHITECTURAL ASSUMPTIONS & EXPERIMENTAL DESIGN DOCUMENTATION
# ==============================================================================
# 1. MODEL GEOMETRY: 2D mid-plane approximation (z = 0) utilizing a 
#    reconstructed radial gravitational acceleration field:
#        a_r(r) = - (V_bar(r)^2 / r) * r_hat
#
# 2. TEST PARTICLES: 200 non-interacting test particles initialized with 
#    logarithmic radial spacing between 1.8 and 42.0 kpc.
#
# 3. EXPERIMENTAL PURPOSE: Test particles are initialized at observed 
#    velocities V_obs(r_i) from SPARC data under BARYONIC GRAVITY ONLY.
#    Because V_obs > V_bar in outer regions, particles are initialized out 
#    of equilibrium to investigate orbital response under baryonic mass deficiency.
# ==============================================================================

# ==============================================================================
# 1. SETUP & DIRECTORIES
# ==============================================================================
BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

CSV_OUTPUT_PATH = RESULTS_DIR / "star_trajectories_baryonic.csv"
SPARC_URL = "https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table2.dat"
VIZIER_URL = (
    "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?"
    "-source=J/AJ/152/157/table2&"
    "-out=Name,Dist,Rad,Vobs,e_Vobs,Vgas,Vdisk,Vbulge,SBdisk,SBbulge&"
    "-out.max=50000"
)
LOCAL_DATA_PATH = BASE_DIR / "table2.dat"

# Physical & Unit Conversion Constants
KM_S_TO_KPC_MYR = 0.00102271  # 1 km/s = 0.00102271 kpc/Myr
UPSILON_DISK = 0.50          # Stellar disk Mass-to-Light ratio
UPSILON_BULGE = 0.70         # Bulge Mass-to-Light ratio

# Simulation Parameters
N_STARS = 200                 # Number of test particles
TIME_TOTAL_MYR = 300.0        # Total simulation time in Myr
N_FRAMES = 301                # 301 keyframes = exactly 1.0 Myr dt per frame
R_MIN_KPC = 1.8               # Inner sampling radius boundary
R_MAX_KPC = 42.0              # Outer sampling radius boundary

# ==============================================================================
# 2. DATA LOADING, INPUT VALIDATION & PCHIP INTERPOLATION
# ==============================================================================
def _download_data() -> str:
    """Download the SPARC table, falling back when the CDS FTP endpoint is blocked."""
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    for url in (SPARC_URL, VIZIER_URL):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, context=ctx, timeout=60) as response:
                text = response.read().decode("utf-8", errors="replace")
            if "NGC3198" in text.upper() and "<html" not in text[:1000].lower():
                return text
            print(f"[INFO] Ignoring invalid response from {url}")
        except (urllib.error.URLError, TimeoutError) as error:
            print(f"[INFO] Download failed from {url}: {error}")
    raise RuntimeError("Unable to download a valid SPARC table from CDS or VizieR.")


def _parse_sparc_rows(text: str) -> list[list[str]]:
    """Parse either the raw CDS table or VizieR's TSV response."""
    rows = []
    for line in text.splitlines():
        parts = line.strip().split()
        if not parts:
            continue
        if parts[0].isdigit() and len(parts) >= 11:
            parts = parts[1:]
        if len(parts) >= 10:
            clean_name = parts[0].replace("-", "").replace("_", "").upper()
            if clean_name == "NGC3198":
                rows.append(parts[:10])
    return rows


if LOCAL_DATA_PATH.exists():
    existing_text = LOCAL_DATA_PATH.read_text(encoding="utf-8", errors="ignore")
    if "<html" in existing_text[:1000].lower() or "NGC3198" not in existing_text.upper():
        print("[INFO] Invalid local table2.dat detected; downloading a fresh copy...")
        LOCAL_DATA_PATH.unlink()

if LOCAL_DATA_PATH.exists():
    table_text = LOCAL_DATA_PATH.read_text(encoding="utf-8", errors="ignore")
else:
    print("[INFO] Downloading SPARC rotation data table...")
    table_text = _download_data()
    LOCAL_DATA_PATH.write_text(table_text, encoding="utf-8")

ngc3198_rows = _parse_sparc_rows(table_text)

if not ngc3198_rows:
    raise ValueError("SPARC data error: No valid NGC3198 entries found in table2.dat.")

# Construct DataFrame
cols = ['Name', 'Dist', 'Rad', 'Vobs', 'eVobs', 'Vgas', 'Vdisk', 'Vbul', 'SBdisk', 'SBbul']
data = pd.DataFrame(ngc3198_rows, columns=cols)

# Convert numerical columns and clean indexing
for col in ['Rad', 'Vobs', 'eVobs', 'Vgas', 'Vdisk', 'Vbul']:
    if col in data.columns:
        data[col] = pd.to_numeric(data[col], errors='coerce')

data = data.dropna(subset=['Rad', 'Vobs', 'eVobs', 'Vgas', 'Vdisk', 'Vbul']).reset_index(drop=True)

# Input Data Validation
raw_columns = ['Rad', 'Vobs', 'eVobs', 'Vgas', 'Vdisk', 'Vbul']
if not np.all(np.isfinite(data[raw_columns].values)):
    raise ValueError("SPARC data validation error: Non-finite values (NaN/Inf) detected in dataset.")

r_data = data['Rad'].to_numpy(dtype=float)
if not np.all(np.diff(r_data) > 0):
    raise ValueError("SPARC data validation error: Radii array 'Rad' is not strictly increasing.")

# Compute total baryonic velocity: V_bar = sqrt(Vgas^2 + Y_d*Vdisk^2 + Y_b*Vbul^2)
v_gas_sq = np.sign(data['Vgas']) * (data['Vgas'] ** 2)
v_disk_sq = UPSILON_DISK * (data['Vdisk'] ** 2)
v_bul_sq = UPSILON_BULGE * (data['Vbul'] ** 2)
v_bar_sq = v_gas_sq + v_disk_sq + v_bul_sq
data['Vbar'] = np.sqrt(np.maximum(0.0, v_bar_sq))

vbar_data = data['Vbar'].values
vobs_data = data['Vobs'].values

r_min_obs = float(r_data[0])
r_max_obs = float(r_data[-1])

# PCHIP Interpolators
pchip_vbar = PchipInterpolator(r_data, vbar_data)
pchip_vobs = PchipInterpolator(r_data, vobs_data)

def get_vbar(r: float) -> float:
    """
    Evaluates V_bar(r) using PCHIP interpolation with explicit boundary domain physics:
    - Inner region (r < r_min_obs): Linear scaling from origin.
    - Observational domain (r_min_obs <= r <= r_max_obs): PCHIP monotonic curve.
    - Outer region (r > r_max_obs): Keplerian falloff V_bar ~ 1 / sqrt(r).
    """
    if r <= 0.0:
        return 0.0
    elif r < r_min_obs:
        return float(pchip_vbar(r_min_obs)) * (r / r_min_obs)
    elif r <= r_max_obs:
        return float(pchip_vbar(r))
    else:
        v_edge = float(pchip_vbar(r_max_obs))
        return v_edge * np.sqrt(r_max_obs / r)

def get_vobs(r: float) -> float:
    """
    Evaluates V_obs(r) using PCHIP interpolation with flat extension outside domain.
    """
    if r <= 0.0:
        return 0.0
    elif r < r_min_obs:
        return float(pchip_vobs(r_min_obs)) * (r / r_min_obs)
    elif r <= r_max_obs:
        return float(pchip_vobs(r))
    else:
        return float(pchip_vobs(r_max_obs))

# ==============================================================================
# 3. GRAVITATIONAL POTENTIAL & ESCAPE VELOCITY ENGINE
# ==============================================================================
# Analytical boundary potential at r_max_obs for a Keplerian outer tail:
# Phi(r) = - (V_bar(r)^2) for r >= r_max_obs, approaching 0 as r -> infinity.
v_edge_kpc_myr = get_vbar(r_max_obs) * KM_S_TO_KPC_MYR
phi_edge = -(v_edge_kpc_myr ** 2)

# High-density numerical grid for interior potential Phi(r) <= r_max_obs
r_pot_grid = np.linspace(0.01, r_max_obs, 2000)
vbar_grid_kpc_myr = np.array([get_vbar(r) for r in r_pot_grid]) * KM_S_TO_KPC_MYR
accel_grid = (vbar_grid_kpc_myr ** 2) / r_pot_grid

dr_pot = np.diff(r_pot_grid)
trapz_areas = 0.5 * (accel_grid[:-1] + accel_grid[1:]) * dr_pot

phi_grid = np.zeros_like(r_pot_grid)
# Integrate inward from r_max_obs: Phi(r) = Phi_edge - integral_r^{r_max_obs} a_r(r') dr'
phi_grid[:-1] = phi_edge - np.cumsum(trapz_areas[::-1])[::-1]
phi_grid[-1] = phi_edge

pchip_phi = PchipInterpolator(r_pot_grid, phi_grid)

def get_phi(r: float) -> float:
    """
    Returns specific gravitational potential Phi(r) in (kpc/Myr)^2.
    Normalized to zero at infinity by incorporating the analytical Keplerian tail.
    Uses quadratic harmonic continuation for r <= 0.01 kpc.
    """
    if r <= 0.01:
        r0 = 0.01
        v0 = get_vbar(r0) * KM_S_TO_KPC_MYR
        phi0 = float(pchip_phi(r0))

        return phi0 + 0.5 * (v0 / r0)**2 * (r**2 - r0**2)
    elif r <= r_max_obs:
        return float(pchip_phi(r))
    else:
        v_kms = get_vbar(r)
        v_kpc_myr = v_kms * KM_S_TO_KPC_MYR
        return -(v_kpc_myr ** 2)

def get_vesc_kms(r: float) -> float:
    """Returns local baryonic escape velocity in km/s: V_esc = sqrt(2 * |Phi(r)|)."""
    phi_val = get_phi(r)
    return np.sqrt(2.0 * abs(phi_val)) / KM_S_TO_KPC_MYR

# ==============================================================================
# 4. 2D MID-PLANE EQUATIONS OF MOTION
# ==============================================================================
def baryonic_equations_of_motion(t, state):
    x, y, vx, vy = state
    r = np.sqrt(x**2 + y**2)
    
    # Smooth core singularity protection
    r_eff = max(r, 1e-4)
        
    vbar_kms = get_vbar(r_eff)
    vbar_kpc_myr = vbar_kms * KM_S_TO_KPC_MYR
    
    # Radial acceleration magnitude: a_r = V_bar^2 / r
    accel_mag = (vbar_kpc_myr ** 2) / r_eff
    
    # Planar acceleration components (-r_hat)
    ax = -accel_mag * (x / r_eff)
    ay = -accel_mag * (y / r_eff)
    
    return [vx, vy, ax, ay]

# ==============================================================================
# 5. INITIAL CONDITIONS, NUMERICAL INTEGRATION & CONSERVATION AUDIT
# ==============================================================================
np.random.seed(42)  # Reproducible placement
t_eval = np.linspace(0.0, TIME_TOTAL_MYR, N_FRAMES)  # Exactly 301 points (dt = 1.0 Myr)
trajectory_records = []

# Diagnostic Tracking
max_L_relative_deviations = []
max_E_relative_deviations = []
max_E_absolute_deviations = []
inner_extrap_sample_count = 0
outer_extrap_sample_count = 0
unbound_particle_count = 0
total_sample_points = N_STARS * N_FRAMES

print("=" * 80)
print("NGC 3198 BARYONIC ORBIT SIMULATOR (2D MID-PLANE APPROXIMATION)")
print("Numerical Safeguards, Energy Scale Normalization & Full Precision Export Active")
print("=" * 80)

# Sample 200 test particles logarithmically across initial radii
radii = np.geomspace(R_MIN_KPC, R_MAX_KPC, N_STARS)
angles = np.random.uniform(0.0, 2.0 * np.pi, N_STARS)

for star_id in range(N_STARS):
    r_i = radii[star_id]
    theta_i = angles[star_id]
    
    # Initial 2D positions in mid-plane (z = 0)
    x0 = r_i * np.cos(theta_i)
    y0 = r_i * np.sin(theta_i)
    
    # Initial observed orbital velocity launch (km/s -> kpc/Myr)
    v_obs_kms = get_vobs(r_i)
    v_obs_kpc_myr = v_obs_kms * KM_S_TO_KPC_MYR
    
    # Pure tangential initial velocity vector
    vx0 = -v_obs_kpc_myr * np.sin(theta_i)
    vy0 =  v_obs_kpc_myr * np.cos(theta_i)
    
    initial_state = [x0, y0, vx0, vy0]
    
    # Initial Specific Angular Momentum L0 and Initial Specific Mechanical Energy E0
    L0 = x0 * vy0 - y0 * vx0
    phi_0 = get_phi(r_i)
    E0 = 0.5 * (v_obs_kpc_myr ** 2) + phi_0
    
    # Energetically unbound particle check (E0 > 0 when Phi is normalized to 0 at infinity)
    if E0 > 0.0:
        unbound_particle_count += 1
    
    # Solve ODE using high-precision Runge-Kutta 45
    sol = solve_ivp(
        baryonic_equations_of_motion,
        t_span=(0.0, TIME_TOTAL_MYR),
        y0=initial_state,
        t_eval=t_eval,
        method='RK45',
        rtol=1e-8,
        atol=1e-10
    )
    
    # Integration Status & Frame Length Verification
    if sol.status != 0:
        raise RuntimeError(f"Integration failed for particle {star_id}: {sol.message}")
    if len(sol.t) != N_FRAMES:
        raise RuntimeError(f"Incomplete integration for particle {star_id}: Got {len(sol.t)} frames, expected {N_FRAMES}.")
    
    # Trajectory arrays
    x_arr = sol.y[0]
    y_arr = sol.y[1]
    vx_arr = sol.y[2]
    vy_arr = sol.y[3]
    
    # Angular Momentum Conservation Check
    L_t = x_arr * vy_arr - y_arr * vx_arr
    L_dev = np.abs(L_t - L0) / np.abs(L0)
    max_L_relative_deviations.append(np.max(L_dev))
    
    # Mechanical Energy Conservation Check with Energy Scale Normalization
    r_t_kpc = np.sqrt(x_arr**2 + y_arr**2)
    v_t_kpc_myr = np.sqrt(vx_arr**2 + vy_arr**2)
    phi_t = np.array([get_phi(r) for r in r_t_kpc])
    E_t = 0.5 * (v_t_kpc_myr**2) + phi_t
    
    energy_abs_error = np.max(np.abs(E_t - E0))
    energy_scale = 0.5 * (v_obs_kpc_myr ** 2) + abs(phi_0)
    energy_relative_error = (
        energy_abs_error / energy_scale
        if energy_scale > 0
        else energy_abs_error
    )
    
    max_E_absolute_deviations.append(energy_abs_error)
    max_E_relative_deviations.append(energy_relative_error)
    
    # Store trajectory keyframes (FULL UNROUNDED NUMERICAL PRECISION)
    for f_idx in range(N_FRAMES):
        x_f = x_arr[f_idx]
        y_f = y_arr[f_idx]
        vx_f = vx_arr[f_idx]
        vy_f = vy_arr[f_idx]
        
        r_f = r_t_kpc[f_idx]
        v_current_kms = v_t_kpc_myr[f_idx] / KM_S_TO_KPC_MYR
        
        # Extrapolation Domain Occupancy Sample Tracking
        if r_f < r_min_obs:
            inner_extrap_sample_count += 1
        elif r_f > r_max_obs:
            outer_extrap_sample_count += 1
            
        trajectory_records.append({
            'frame': f_idx + 1,
            'time_myr': t_eval[f_idx],
            'star_id': star_id,
            'initial_r_kpc': r_i,
            'v_obs_kms': v_obs_kms,
            'current_v_kms': v_current_kms,
            'x_kpc': x_f,
            'y_kpc': y_f,
            'z_kpc': 0.0
        })

# ==============================================================================
# 6. DIAGNOSTICS REPORTING & CSV EXPORT
# ==============================================================================
# Export full precision dataset to CSV
df_trajectories = pd.DataFrame(trajectory_records)
df_trajectories.to_csv(CSV_OUTPUT_PATH, index=False)

# Compute Conservation & Domain Occupancy Statistics
avg_max_L_dev = np.mean(max_L_relative_deviations)
worst_L_dev = np.max(max_L_relative_deviations)
avg_max_E_rel_dev = np.mean(max_E_relative_deviations)
worst_E_rel_dev = np.max(max_E_relative_deviations)

pct_inner_extrap = (inner_extrap_sample_count / total_sample_points) * 100.0
pct_outer_extrap = (outer_extrap_sample_count / total_sample_points) * 100.0

vesc_initial_kms = [get_vesc_kms(r) for r in radii]

print(f"[SUCCESS] Exported {len(df_trajectories)} unrounded frame keyframes to:")
print(f"         {CSV_OUTPUT_PATH}")
print("-" * 80)
print("NUMERICAL CONSERVATION AUDIT:")
print(f" - Avg Max Relative Angular Momentum Deviation (|dL/L0|): {avg_max_L_dev:.4e}")
print(f" - Worst Particle Relative Angular Momentum Deviation:    {worst_L_dev:.4e}")
print(f" - Avg Max Scaled Relative Energy Error (dE / E_scale):    {avg_max_E_rel_dev:.4e}")
print(f" - Worst Particle Scaled Relative Energy Error:            {worst_E_rel_dev:.4e}")
print("-" * 80)
print("ROTATION CURVE DOMAIN OCCUPANCY TRACKING:")
print(f" - Main SPARC Observed Domain ({r_min_obs:.2f} to {r_max_obs:.2f} kpc): {100.0 - (pct_inner_extrap + pct_outer_extrap):.2f}% of samples")
print(f" - Inner Linear Extrapolation Domain (< {r_min_obs:.2f} kpc):          {pct_inner_extrap:.2f}% of samples")
print(f" - Outer Keplerian Falloff Domain (> {r_max_obs:.2f} kpc):           {pct_outer_extrap:.2f}% of samples")
print("-" * 80)
print("KINEMATIC OUTCOME SUMMARY:")
print(f" - Initial Baryonic Escape Velocity Range:                  {min(vesc_initial_kms):.1f} to {max(vesc_initial_kms):.1f} km/s")
print(f" - Energetically unbound particles:                        {unbound_particle_count} / {N_STARS} ({(unbound_particle_count/N_STARS)*100.0:.1f}%)")
print("=" * 80)