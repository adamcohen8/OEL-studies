"""Verify the final cover, global chart limits, witness membership, and detector numbers."""
from collections import Counter
from hashlib import sha256
from math import atan, atan2, asin, degrees, radians, sqrt, tan
from audit_partition import CachedDual, exact, install, sc, GG, gs, limits, D, np, json
from outer_test import path, n


def read(name):
    return json.loads((D / (name + '.json')).read_text())


def collect(record, box):
    if 'children' not in record:
        assert np.array_equal(record['box'], box)
        yield box, record['result']
        return
    axis, cut = record['axis'], record['cut']
    assert box[axis, 0] < cut < box[axis, 1] and len(record['children']) == 2
    for side, child in enumerate(record['children']):
        piece = box.copy()
        piece[axis, 1 - side] = cut
        yield from collect(child, piece)


def main():
    install()
    # Assess the root limits against the incoming optical polytope, with no chart constraints.
    solver = CachedDual()
    root_slacks = []
    for axis in range(6):
        for sign in (-1, 1):
            goal = sign * GG[axis]
            assert solver(sc.fullM, sc.fullb, sc.lb, sc.ub, goal) is not None
            upper = exact(solver, sc.fullM, sc.fullb, goal)
            assert upper is not None
            root_slacks.append(float(sign * limits[axis, int(sign > 0)] - upper * gs[axis]))
    assert min(root_slacks) >= 0, root_slacks
    boxes = read('audited-boxes')
    digest = sha256((D / 'audited-boxes.json').read_bytes()).hexdigest()
    groups = []
    saved = set()
    for shard in range(2):
        packet = read(f'final-shard-{shard}')
        assert packet['boxes_sha256'] == digest and packet['shards'] == 2
        assert sorted(packet['records']) == sorted(boxes)[shard::2]
        for key, record in packet['records'].items():
            assert key not in saved
            saved.add(key)
            groups.extend(collect(record, np.array(boxes[key])))
    leaves = {r['id'] for r in read('optical-cover')['leaves']}
    assert saved == leaves == set(boxes)
    routes = read('partition-routes')
    assert len(routes) == len(leaves) - 32
    def visit(key):
        if key not in routes:
            assert key in leaves
            return {key}
        assert key not in leaves and 0 <= routes[key]['axis'] < 6
        return visit(key + '.0') | visit(key + '.1')
    assert set.union(*(visit(str(i)) for i in range(32))) == leaves
    counts = Counter(r['status'] for _, r in groups)
    assert set(counts) <= {'contained', 'infeasible', 'collision', 'all intersect surface'}
    contained = [(b, r) for b, r in groups if r['status'] == 'contained']
    max_excess = max(r['max_outside_plane_km'] for _, r in contained)
    assert max_excess < 0
    chart = read('chart')
    nominal = read('nominal-certificate')['2']
    witness_counts = []
    for witness in nominal['accepted']:
        last = path(np.array(witness['state']))[6*n-6:6*n]
        z = np.array(chart['transform']) @ (last - chart['origin'])
        witness_counts.append(sum(bool(np.all(b[:, 0] <= z) and np.all(z <= b[:, 1])) for b, _ in contained))
    assert min(witness_counts) >= 1
    fit = read('outer-edges')['fit']
    c = fit['center']
    rolls = read('all-rolls')
    group = rolls['2.15']
    assert read('interval-roll')['2.15']['success'] and read('nominal-lower-interval')['success']
    cap = group['max_training_error_as']
    assert cap < 2.148
    circle_upper = degrees(2*atan(sqrt(2)*tan(radians(2.2))))
    assert group['all_pair_distances_upper_deg'] < circle_upper
    result = dict(epoch_utc='2024-08-20T06:00:00Z', observatory_code='703',
                  observation_count=25, nominal_cap_arcsec=2.0,
                  square_side_bounds_deg=[rolls['2']['required_side_lower_deg'], 4.4],
                  circle_diameter_bounds_deg=[nominal['minimum_pair_deg'], circle_upper],
                  boresight=c, ra_deg=degrees(atan2(c[1], c[0])) % 360, dec_deg=degrees(asin(c[2])),
                  roll_east_toward_north_deg=fit['roll_deg'],
                  incoming_slabs=32, internal_splits=len(routes), audited_leaves=len(leaves),
                  unconstrained_chart_projection_slack=root_slacks,
                  chart_supports_rationally_replayed=12, partition_supports_rationally_replayed=read('partition-audit')['exact_support_rechecks'],
                  audited_boxes_sha256=digest, final_classification=counts, unresolved=0,
                  worst_contained_face_upper_km=max_excess, contained_box_counts_for_nominal_witnesses=witness_counts,
                  sufficient_error_cap_arcsec=2.0, excluded_error_cap_arcsec_rounded_up=2.148,
                  three_orbit_maximum_residual_upper_arcsec=cap, three_orbit_indices=group['indices'],
                  three_orbit_max_pair_upper_deg=group['all_pair_distances_upper_deg'],
                  independent_integrator_max_training_residual_arcsec={k:max(r['max_fit_as'] for r in v.values()) for k,v in read('independent').items()})
    (D / 'RESULTS.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
