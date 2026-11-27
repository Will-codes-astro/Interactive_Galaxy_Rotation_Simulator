"""1-sigma uncertainty of every value the site's best fits tune, from the profile chi-squared.

Will's choice (27 Nov 2026): hold one fitted value at a trial value and refit all the others; its 1-sigma range is
where chi-squared has risen by 1 above the minimum, with the error bars scaled up by sqrt(chi2_r) when the fit is poor
(chi2_r > 1), so the rise is max(1, chi2_r). A range that reaches a realistic limit stops there. A halo range end
where the refitted halo sits on its concentration-for-mass limit or the scale radius slider's end is recorded too
(uncertainty_halo_limit): without that limit the range could reach further.

Each refit uses SciPy's global optimiser (differential evolution). Distance, tilt and M/L trade off along long curved
valleys in chi-squared, and the generator's local search, moving one value at a time, can stop short in them: on
27 Nov 2026 it missed fits that differential evolution found. A global search takes minutes, so it runs here, not in
the generator, which reads results/fit_uncertainties.json and stops if the file no longer matches its fits.

Uses the generator's own fits, models and chi-squared (generate_interactive_html.py) and SPARC via sparc_data.py.
Run from the project root: python analysis/fit_uncertainties.py
"""

import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor

from scipy.optimize import differential_evolution

from sparc_data import ROOT, load_galaxy_info

sys.path.insert(0, str(ROOT))
import generate_interactive_html as site  # noqa: E402  (needs ROOT on the path)

OUTPUT_PATH = ROOT / 'results' / 'fit_uncertainties.json'
FIT_KINDS = ('nfw_standard', 'nfw', 'mond')
SETTING_NAMES = site.SETTING_NAMES
# Each range end is bisected until it is known to PROFILE_PRECISION of its distance from the best value (the site shows
# the ranges to 2 significant figures), in at most PROFILE_MAX_BISECTIONS halvings. A fixed number of halvings of the
# gap to the limit is not enough: a halo mass range can be 50 times narrower than that gap (found 27 Nov 2026).
PROFILE_PRECISION = 0.01
PROFILE_MAX_BISECTIONS = 40
# Differential evolution: population per value, convergence tolerance, generation cap and a fixed seed (repeatable).
DE_POPSIZE = 15
DE_TOLERANCE = 1e-8
DE_MAX_GENERATIONS = 1000
DE_SEED = 3198
# Fits outside the realistic limits cost this instead of infinity (the optimiser's polish step needs finite values).
OUT_OF_LIMITS_COST = 1e12
# A refitted halo within this fraction of a limit's span from the limit counts as sitting on it. Ends that sat on a
# limit were within 0.001 dex of the concentration band on 27 Nov 2026; the nearest that did not was 0.013 dex away.
HALO_LIMIT_FRACTION = 0.01


def fit_problem(kind, rotation, info):
	"""The fit as an optimisation problem: value names, limits, chi-squared, each value in the record's units, and
	which of the halo's own limits a set of values sits on (None, or always None for MOND).

	The NFW halo is described by log10 M200 and log10 c, the values the Scorecard shows, and c must stay within
	the generator's realistic limit of the Dutton & Maccio (2014) concentration for its mass. That limit and the
	scale radius slider's range depend on M200 and c together, so they are not edges of either value's own range.
	"""
	settings = site._standard_settings(info) if kind == 'nfw_standard' else site._realistic_settings(rotation, info)
	plain = [lambda x: x] * len(SETTING_NAMES)
	if kind == 'mond':
		bounds = [tuple(a / 1e-10 for a in site.MOND_A0_RANGE_MS2), *settings]

		def cost(p):
			return site._chi_squared(rotation, info, p[1:], lambda v, r: site._mond_velocity(v, r, p[0] * 1e-10))

		return ('a0_ms2', *SETTING_NAMES), bounds, cost, [lambda x: x * 1e-10, *plain], lambda p: None

	log_c_limit = site.FIT_LIMIT_SIGMA * site.HALO_CONCENTRATION_SCATTER_DEX
	log_c = lambda m200: math.log10(site._halo_concentration(m200))
	mass_lo, mass_hi = site.HALO_MASS_RANGE_MSUN
	bounds = [(math.log10(mass_lo), math.log10(mass_hi)), (log_c(mass_hi) - log_c_limit, log_c(mass_lo) + log_c_limit), *settings]

	def cost(p):
		m200 = 10 ** p[0]
		scale_radius = site._r200_kpc(m200) / 10 ** p[1]
		if abs(p[1] - log_c(m200)) > log_c_limit or not site.HALO_SCALE_RANGE_KPC[0] <= scale_radius <= site.HALO_SCALE_RANGE_KPC[1]:
			return OUT_OF_LIMITS_COST
		halo = lambda v, r: math.hypot(v, site._nfw_velocity(r, m200, scale_radius))
		return site._chi_squared(rotation, info, p[2:], halo)

	def halo_limit(p):
		offset = abs(p[1] - log_c(10 ** p[0]))
		if offset >= log_c_limit - HALO_LIMIT_FRACTION * 2 * log_c_limit:
			return 'concentration_band'
		lo, hi = site.HALO_SCALE_RANGE_KPC
		scale_radius = site._r200_kpc(10 ** p[0]) / 10 ** p[1]
		if min(scale_radius - lo, hi - scale_radius) <= HALO_LIMIT_FRACTION * (hi - lo):
			return 'scale_radius'
		return None

	unlog = lambda x: 10 ** x
	return ('m200_msun', 'concentration', *SETTING_NAMES), bounds, cost, [unlog, unlog, *plain], halo_limit


def inside(x, lo, hi):
	"""x moved just inside [lo, hi]: differential evolution rejects a start that its own rescaling puts past a limit."""
	margin = 1e-9 * (hi - lo)
	return min(max(x, lo + margin), hi - margin)


def minimise(cost, bounds, held=None, start=None, stop_at=-math.inf):
	"""Lowest chi-squared over every free value, with held = (index, value) fixed. Stops early once at or below stop_at.

	Returns the chi-squared and the full list of values.
	"""
	free = [i for i, (lo, hi) in enumerate(bounds) if hi > lo and (held is None or i != held[0])]
	base = [lo for lo, _ in bounds]
	if held is not None:
		base[held[0]] = held[1]

	def full(q):
		p = base.copy()
		for i, v in zip(free, q):
			p[i] = v
		return p

	if not free:
		return cost(base), base
	result = differential_evolution(
		lambda q: cost(full(q)), [bounds[i] for i in free], popsize=DE_POPSIZE, tol=DE_TOLERANCE,
		maxiter=DE_MAX_GENERATIONS, seed=DE_SEED,
		x0=None if start is None else [inside(start[i], *bounds[i]) for i in free],
		callback=lambda xk, convergence: cost(full(xk)) <= stop_at,
	)
	return result.fun, full(result.x)


def profile_ranges(job):
	"""One fit's ranges: (galaxy id, fit kind, record)."""
	galaxy_id, kind, rotation, info = job
	if kind == 'mond':
		fit = site.fit_mond(rotation, info)
	else:
		fit = site.fit_nfw(rotation, info, site._standard_settings(info) if kind == 'nfw_standard' else None)
	names, bounds, cost, to_record, halo_limit = fit_problem(kind, rotation, info)
	chi_min, best = minimise(cost, bounds)
	free = sum(hi > lo for lo, hi in bounds)
	rise = max(1.0, chi_min / (len(rotation['radius_kpc']) - free))
	target = chi_min + rise

	def within(d, x, anchor):
		chi, values = minimise(cost, bounds, (d, x), anchor, target)
		return chi <= target, values

	ranges, cut, halo_limited = {}, {}, {}
	for d, name in enumerate(names):
		lo, hi = bounds[d]
		if hi <= lo:
			continue
		ends = []
		for side, edge in (('low', lo), ('high', hi)):
			if within(d, edge, best)[0]:
				ends.append(edge)
				cut[name] = 'both' if name in cut else side
				continue
			inside, outside, anchor = best[d], edge, best
			for _ in range(PROFILE_MAX_BISECTIONS):
				if abs(outside - inside) <= PROFILE_PRECISION * abs(inside - best[d]):
					break
				middle = (inside + outside) / 2.0
				accepted, values = within(d, middle, anchor)
				if accepted:
					inside, anchor = middle, values
				else:
					outside = middle
			ends.append((inside + outside) / 2.0)
			# An end can also be set by a halo limit that is not an edge of this value's range: refit fully (the
			# bisection stops each refit once it is good enough) and see whether the best halo there sits on one.
			if name in ('m200_msun', 'concentration'):
				limit = halo_limit(minimise(cost, bounds, (d, inside), anchor)[1])
				if limit:
					halo_limited.setdefault(name, {})[side] = limit
		ranges[name] = [float(f'{to_record[d](x):.4g}') for x in ends]
	return galaxy_id, kind, {
		'fit_chi2': fit['chi2'],
		'chi2_min': round(chi_min, 3),
		'chi2_rise': round(rise, 3),
		'uncertainty': ranges,
		'uncertainty_cut': cut,
		'uncertainty_halo_limit': halo_limited,
	}


def main():
	info = load_galaxy_info()
	jobs = []
	for preset in site.GALAXY_PRESETS:
		rotation = site.load_sparc_galaxy(preset.sparc_id)['rotation']
		jobs += [(preset.sparc_id, kind, rotation, info[preset.sparc_id.upper()]) for kind in FIT_KINDS]
	results = {}
	with ProcessPoolExecutor() as pool:
		for galaxy_id, kind, record in pool.map(profile_ranges, jobs):
			results.setdefault(galaxy_id, {})[kind] = record
			print(f'{galaxy_id} {kind}: chi2 min {record["chi2_min"]}, rise {record["chi2_rise"]}, {record["uncertainty"]}', flush=True)
	OUTPUT_PATH.write_text(json.dumps({
		'method': 'profile chi-squared, rise max(1, chi2_r) above the minimum; differential evolution refits; '
		          'ranges stop at the realistic limits (see analysis/fit_uncertainties.py)',
		'galaxies': results,
	}, indent=1) + '\n', encoding='utf-8')
	print(f'Wrote {OUTPUT_PATH}')


if __name__ == '__main__':
	main()
