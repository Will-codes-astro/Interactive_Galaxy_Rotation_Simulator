"""Baryonic-only rotation curve, residuals and reduced chi-squared for SPARC galaxies.

Reads data/sparc/table2.dat (via sparc_data.py) and writes plots/<galaxy>_baryonic_baseline.png/.pdf
and results/rotation_curve_summary.csv/.txt. Run from the project root: python analysis/baryonic_baseline.py
"""

import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sparc_data import ROOT, load_table

BASE_DIR = ROOT
PLOTS_DIR = BASE_DIR / 'plots'
RESULTS_DIR = BASE_DIR / 'results'
PLOTS_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True)

# Add future SPARC galaxy targets to this list.
GALAXIES_TO_ANALYSE = ['NGC3198']


def _normalise_galaxy_name(name: object) -> str:
	"""Canonicalise galaxy names for case-insensitive lookups."""
	return re.sub(r'[^a-z0-9]+', '', str(name).lower())


def load_sparc_data() -> pd.DataFrame:
	"""The validated local SPARC table as a DataFrame (column names as in SPARC)."""
	rows = load_table()
	return pd.DataFrame(
		{
			'Name': [row.name for row in rows],
			'Dist': [row.dist_mpc for row in rows],
			'Rad': [row.rad_kpc for row in rows],
			'Vobs': [row.v_obs_kms for row in rows],
			'e_Vobs': [row.e_v_obs_kms for row in rows],
			'Vgas': [row.v_gas_kms for row in rows],
			'Vdisk': [row.v_disk_kms for row in rows],
			'Vbul': [row.v_bulge_kms for row in rows],
		}
	)


def analyse_galaxy(
	data: pd.DataFrame,
	galaxy_name: str,
	ups_disk: float = 0.5,
	ups_bul: float = 0.7,
) -> dict[str, object]:
	"""Calculate, plot, and export one galaxy's baryonic baseline."""
	normalized_target = _normalise_galaxy_name(galaxy_name)
	galaxy = data[data['Name'].map(_normalise_galaxy_name) == normalized_target].copy()
	if galaxy.empty:
		raise ValueError(f'No SPARC data found for {galaxy_name}.')

	galaxy['Vbar'] = np.sqrt(
		ups_disk * galaxy['Vdisk'] ** 2
		+ ups_bul * galaxy['Vbul'] ** 2
		+ galaxy['Vgas'] ** 2
	)
	residuals = galaxy['Vobs'] - galaxy['Vbar']
	valid_errors = galaxy['e_Vobs'].gt(0) & galaxy['Vobs'].notna() & galaxy['Vbar'].notna()
	if not valid_errors.any():
		raise ValueError(f'No positive observational errors found for {galaxy_name}.')
	chi_squared = np.sum(
		(
			residuals[valid_errors]
			/ galaxy.loc[valid_errors, 'e_Vobs']
		) ** 2
	)
	degrees_of_freedom = int(valid_errors.sum())
	if degrees_of_freedom == 0:
		raise ValueError(f'No valid measurements available for {galaxy_name}.')
	reduced_chi_squared = chi_squared / degrees_of_freedom
	outer_point = galaxy.loc[galaxy['Rad'].idxmax()]
	max_baryonic_velocity = galaxy['Vbar'].max()
	max_observed_velocity = galaxy['Vobs'].max()
	outer_velocity_deficit = outer_point['Vobs'] - outer_point['Vbar']

	filename = re.sub(r'[^a-z0-9]+', '_', galaxy_name.lower()).strip('_')
	png_path = PLOTS_DIR / f'{filename}_baryonic_baseline.png'
	pdf_path = PLOTS_DIR / f'{filename}_baryonic_baseline.pdf'

	figure, (axis, residual_axis) = plt.subplots(
		2,
		1,
		figsize=(9, 9),
		sharex=True,
		gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.08},
	)
	axis.errorbar(
		galaxy['Rad'],
		galaxy['Vobs'],
		yerr=galaxy['e_Vobs'],
		fmt='o',
		color='black',
		ecolor='gray',
		capsize=3,
		label='Observed Velocity ($V_{obs}$)',
	)
	axis.plot(
		galaxy['Rad'],
		galaxy['Vgas'],
		':',
		color='teal',
		label='Gas Component ($V_{gas}$)',
	)
	axis.plot(
		galaxy['Rad'],
		np.sqrt(ups_disk) * galaxy['Vdisk'],
		':',
		color='crimson',
		label='Disk Component ($V_{disk}$)',
	)
	axis.plot(
		galaxy['Rad'],
		galaxy['Vbar'],
		'--',
		color='blue',
		linewidth=2,
		label='Total Baryonic Baseline ($V_{bar}$)',
	)
	axis.set_ylabel('Rotation Velocity $V$ (km/s)', fontsize=12)
	axis.set_title(
		f'{galaxy_name}: Observed Rotation vs. Baryonic Expectations',
		fontsize=14,
		pad=12,
	)
	axis.legend(frameon=True, facecolor='white', framealpha=0.9)
	axis.grid(True, linestyle='--', alpha=0.5)
	axis.text(
		galaxy['Rad'].max(),
		100,
		(
			f'Mass-to-light ratios (solar units):\n'
			f'  disk = {ups_disk:.2f}, bulge = {ups_bul:.2f}\n'
			f'Max baryonic velocity = {max_baryonic_velocity:.2f} km/s\n'
			f'Max observed velocity = {max_observed_velocity:.2f} km/s\n'
			f'Outer velocity deficit = {outer_velocity_deficit:.2f} km/s\n'
			f'Reduced chi-squared = {reduced_chi_squared:.2f}\n'
			f'(positive deficit means baryons fall short)'
		),
		transform=axis.transData,
		ha='right',
		va='bottom',
		fontsize=8.5,
		bbox={
			'boxstyle': 'round,pad=0.5',
			'facecolor': 'white',
			'edgecolor': 'gray',
			'alpha': 0.78,
		},
	)
	residual_axis.errorbar(
		galaxy['Rad'],
		residuals,
		yerr=galaxy['e_Vobs'],
		fmt='o',
		color='darkorange',
		ecolor='gray',
		capsize=3,
		label='$V_{obs} - V_{bar}$',
	)
	residual_axis.axhline(0, color='black', linewidth=1, linestyle='--')
	residual_axis.set_xlabel('Galactocentric Radius $R$ (kpc)', fontsize=12)
	residual_axis.set_ylabel('Residual\n(km/s)', fontsize=10)
	residual_axis.set_title('Residuals: Observed minus Baryonic Velocity', fontsize=11)
	residual_axis.grid(True, linestyle='--', alpha=0.5)
	residual_axis.legend(frameon=True, facecolor='white', framealpha=0.9)
	figure.align_ylabels((axis, residual_axis))
	figure.savefig(png_path, dpi=300, bbox_inches='tight')
	figure.savefig(pdf_path, format='pdf', bbox_inches='tight')
	plt.close(figure)

	print(
		f'{galaxy_name}: reduced Chi-squared = {reduced_chi_squared:.2f}, '
		f'outer velocity deficit = {outer_velocity_deficit:.2f} km/s'
	)
	print(f'Exported {png_path.name} and {pdf_path.name}')
	return {
		'Galaxy': galaxy_name,
		'Distance_Mpc': galaxy['Dist'].iloc[0],
		'Data_Points': len(galaxy),
		'Ups_disk': ups_disk,
		'Ups_bul': ups_bul,
		'Degrees_of_Freedom': degrees_of_freedom,
		'Reduced_Chi2': reduced_chi_squared,
		'Max_Baryonic_Velocity_km_s': max_baryonic_velocity,
		'Max_Observed_Velocity_km_s': max_observed_velocity,
		'Outer_Radius_kpc': outer_point['Rad'],
		'Outer_Observed_Velocity_km_s': outer_point['Vobs'],
		'Outer_Baryonic_Velocity_km_s': outer_point['Vbar'],
		'Outer_Velocity_Deficit_km_s': outer_velocity_deficit,
		'Plot_PNG': str(png_path.relative_to(BASE_DIR)),
		'Plot_PDF': str(pdf_path.relative_to(BASE_DIR)),
	}


def main() -> None:
	data = load_sparc_data()
	fitting_results = [
		analyse_galaxy(data, galaxy_name)
		for galaxy_name in GALAXIES_TO_ANALYSE
	]
	results_df = pd.DataFrame(fitting_results)
	numeric_columns = results_df.select_dtypes(include='number').columns
	results_df[numeric_columns] = results_df[numeric_columns].astype(float)
	results_path = RESULTS_DIR / 'rotation_curve_summary.csv'
	results_df.to_csv(results_path, index=False, float_format='%.2f')
	print(f'Results summary exported to {results_path}')

	plain_english_lines = [
		'Plain-English Rotation-Curve Summary',
		'====================================',
		'',
		'The script reads the local SPARC table (data/sparc/table2.dat). The table is',
		'read as whitespace-separated columns with galaxy name, distance, radius,',
		'observed rotation velocity, its observational error, and the gas, disk,',
		'and bulge velocity components. Only the requested galaxy names are then',
		'selected for analysis.',
		'',
		'For each galaxy, the baryonic velocity is calculated by adding the squared',
		'gas, disk, and bulge contributions in quadrature. The disk and bulge',
		'contributions are scaled by their mass-to-light ratios before this sum:',
		'Vbar = sqrt(Vgas^2 + Ups_disk*Vdisk^2 + Ups_bul*Vbul^2).',
		'The adopted mass-to-light ratios are fixed assumptions, not fitted values.',
		'',
		'Reduced chi-squared compares the observed and baryonic velocities after',
		'dividing each difference by its quoted observational error. It is the',
		'sum of those squared, error-weighted differences divided by the number',
		'of valid measurements. A value near 1 indicates agreement at the level',
		'of the reported errors; a much larger value indicates that the baryonic',
		'profile fits the data poorly relative to those errors.',
		'',
	]
	for result in fitting_results:
		deficit = float(result['Outer_Velocity_Deficit_km_s'])
		if deficit > 0:
			outer_explanation = (
				'The positive deficit means the observed outer rotation is faster '
				'than the baryonic prediction.'
			)
		else:
			outer_explanation = (
				'The negative deficit means the baryonic prediction is at least as '
				'fast as the observed outer rotation.'
			)
		plain_english_lines.extend([
			f"{result['Galaxy']}",
			'-' * len(str(result['Galaxy'])),
			(
				f"The maximum baryonic velocity is "
				f"{float(result['Max_Baryonic_Velocity_km_s']):.2f} km/s, while "
				f"the maximum observed velocity is "
				f"{float(result['Max_Observed_Velocity_km_s']):.2f} km/s."
			),
			(
				f"At the outermost measured radius of "
				f"{float(result['Outer_Radius_kpc']):.2f} kpc, the observed velocity "
				f"minus the baryonic velocity is {deficit:.2f} km/s. "
				f"{outer_explanation}"
			),
			(
				f"The reduced chi-squared is "
				f"{float(result['Reduced_Chi2']):.2f}, using "
				f"{float(result['Degrees_of_Freedom']):.2f} valid measurements."
			),
			'',
		])
	summary_text_path = RESULTS_DIR / 'rotation_curve_summary.txt'
	summary_text_path.write_text('\n'.join(plain_english_lines), encoding='utf-8')
	print(f'Plain-English summary exported to {summary_text_path}')


if __name__ == '__main__':
	main()

