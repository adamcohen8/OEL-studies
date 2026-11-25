import json
from pathlib import Path

r = json.loads(Path('results.json').read_text(encoding='utf-8'))
rows = []

def compact(x):
    if isinstance(x, float):
        return float(format(x, '.8g'))
    if isinstance(x, dict):
        return {k: compact(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [compact(v) for v in x]
    return x

def add(label, **values):
    rows.append([label, json.dumps(compact(values), separators=(',', ':'), allow_nan=False)])

def grid_id(g):
    return f"{g['r0_au']:g} AU beta={g['beta']:g}"

add('Environment', **r['environment'])
f = r['quadrature_fit']
add('Cutoff-corrected quadrature fit', analytic_K_yr=r['analytic_K_years'], **f)
grid = r['grid']
for key in ('quad_coarse_vs_fine_relative', 'quad_fine_vs_secular_relative'):
    g = max(grid, key=lambda g: abs(g[key]))
    add('Grid max absolute ' + key, case=grid_id(g), signed_rel=g[key], secular_yr=g['secular_years'], quadrature_yr=[s['time_years'] for s in g['quadrature']])

def drift_error(node):
    return node['relative_error_vs_minus_2k_over_r']

def phase_span(node):
    rates = node['phase_da_dt_km_s']
    return (max(rates) - min(rates)) / abs(node['mean_da_dt_km_s'])

nodes = [(g, s, n) for g in grid for s in g['quadrature'] for n in s['nodes']]
for label, metric in [('mean drift relative error', drift_error), ('relative phase span', phase_span)]:
    g, s, n = max(nodes, key=lambda t: abs(metric(t[2])))
    add('Quadrature node max absolute ' + label, case=grid_id(g), order=s['order'], phases=s['phases'], radius_km=n['radius_km'], value=metric(n), mean_dadt_km_s=n['mean_da_dt_km_s'])

ratio = 'measured_to_secular_decay_ratio'
change = 'negative_change_in_a_squared_over_r0_squared'
def tight(g):
    return g['cartesian_windows'][-1]

for label, g in [
    ('weakest signal', min(grid, key=lambda g: tight(g)['expected_negative_change'])),
    ('largest signal', max(grid, key=lambda g: tight(g)['expected_negative_change'])),
    ('largest secular discrepancy', max(grid, key=lambda g: abs(tight(g)[ratio] - 1))),
    ('largest refinement shift', max(grid, key=lambda g: abs(g['cartesian_windows'][0][ratio] - tight(g)[ratio])))
]:
    c, t = g['cartesian_windows']
    add('Short window: ' + label, case=grid_id(g), duration_yr=t['duration_years'], expected_delta_a2=t['expected_negative_change'], measured_delta_a2=t[change], fine_ratio=t[ratio], coarse_minus_fine_ratio=c[ratio]-t[ratio])

for s in r['zero_drag_control']:
    add('Zero drag (1 AU, 16 initial periods)', rtol=s['rtol'], duration_yr=s['duration_years'], signed_delta_a2=s[change], torque_resid=s['max_scaled_angular_momentum_identity_residual'], max_step_rad=s['max_step_angle_rad'])

for c in r['full_cases']:
    tag = grid_id(c) + ' ' + c['launch']
    d, p = c['direct'][-1], c['polar'][-1]
    x = c['comparisons'][-1]
    add('Surface: ' + tag, q=x['q'], q_squared=c['omitted_fraction_if_surface_time_is_approximated_by_zero_cutoff_law'], direct_yr=d['events'][-1]['time_years'], polar_yr=p['events'][-1]['time_years'], secular_yr=x['circular_secular_time_years'], D_vs_sec=x['direct_vs_secular_relative'], D_minus_sec_yr=x['direct_minus_secular_years'], D_vs_P=x['direct_vs_polar_relative'], D_refine=x['direct_refinement_relative'], P_refine=x['polar_refinement_relative'])
    add('Diagnostics: ' + tag, lambda0=c['lambda0'], lambda_surface_circ=c['lambda_circular_at_solar_surface'], max_v_over_c=d['max_sampled_speed_over_c'], torque_resid=d['max_scaled_angular_momentum_identity_residual'], D_max_step_rad=d['max_step_angle_rad'], D_P_steps=[d['accepted_steps'],p['accepted_steps']], D_P_surface_vr_km_s=[s['events'][-1]['radial_velocity_km_s'] for s in (d,p)])
    inner = c['comparisons'][:-1]
    def peak(key):
        row = max(inner, key=lambda row: abs(row[key]))
        return [row['q'], row[key]]
    add('Interior cuts: ' + tag, q=[x['q'] for x in inner], D_yr=[e['time_years'] for e in d['events'][:-1]], D_vs_sec=[x['direct_vs_secular_relative'] for x in inner], peak_D_vs_P_q_rel=peak('direct_vs_polar_relative'), peak_D_refine_q_rel=peak('direct_refinement_relative'), peak_P_refine_q_rel=peak('polar_refinement_relative'))

repeats = [(grid_id(c), s['method'], s['rtol'], e['q'], e['crossing_count']) for c in r['full_cases'] for s in c['direct'] + c['polar'] for e in s['events'] if e['crossing_count'] > 1]
if repeats:
    add('Repeated inward crossings', count=len(repeats), largest_count=max(repeats, key=lambda x: x[-1]))

Path('summary.json').write_text(json.dumps(rows, separators=(',', ':')) + '\n', encoding='utf-8')
