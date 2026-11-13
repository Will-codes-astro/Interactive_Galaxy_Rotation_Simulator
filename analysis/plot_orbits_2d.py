import pandas as pd
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from pathlib import Path

# Load trajectory data
data_path = Path(__file__).resolve().parent.parent / "results" / "star_trajectories_baryonic.csv"
if not data_path.exists():
    raise FileNotFoundError("Could not find star_trajectories_baryonic.csv in results folder.")

df = pd.read_csv(data_path)

# Set up the plot window
fig, ax = plt.subplots(figsize=(8, 8), facecolor='black')
ax.set_facecolor('black')

# Set plot bounds (kpc)
max_r = df['initial_r_kpc'].max() * 1.1
ax.set_xlim(-max_r, max_r)
ax.set_ylim(-max_r, max_r)
ax.set_aspect('equal')

# Title & labels styling
ax.tick_params(colors='white')
for spine in ax.spines.values():
    spine.set_color('#444444')

ax.set_xlabel("x (kpc)", color='white')
ax.set_ylabel("y (kpc)", color='white')

# Scatter plot for stars
stars = ax.scatter([], [], c='#00e5ff', s=12, alpha=0.8, edgecolors='none')
time_text = ax.text(0.03, 0.95, '', transform=ax.transAxes, color='white', fontsize=12)

# Frame update function
frames = sorted(df['frame'].unique())

def update(frame_num):
    frame_data = df[df['frame'] == frame_num]
    stars.set_offsets(frame_data[['x_kpc', 'y_kpc']].values)
    time_val = frame_data['time_myr'].iloc[0]
    time_text.set_text(f"NGC 3198 (Baryonic Only)\nTime: {time_val:.1f} Myr")
    return stars, time_text

plt.title("NGC 3198 Test Particle Orbits", color='white', pad=15)
plt.tight_layout()

animation_path = Path(__file__).resolve().parent.parent / "plots" / "ngc3198_baryonic_orbits_2d.gif"
animation_path.parent.mkdir(exist_ok=True)
anim = FuncAnimation(fig, update, frames=frames, interval=1000 / 30, blit=True)
anim.save(animation_path, writer=PillowWriter(fps=30))
plt.close(fig)
print(f"[SUCCESS] Saved 2D animation to: {animation_path}")