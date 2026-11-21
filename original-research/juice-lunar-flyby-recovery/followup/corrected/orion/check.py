"""Independent replays of the reported witnesses and review-store truth."""

import json
import numpy as np
import orion as o


from sim.dynamics.orbit.propagator import OrbitPropagator


def angle_visibility(u, observer, moon):
    to_moon = moon - observer
    # Angular clearance relative to the lunar limb (conservative even if Moon behind target).
    return o.distance(u, o.unit(to_moon)) / 3600 - np.rad2deg(np.arcsin(1737.4 / np.linalg.norm(to_moon, axis=-1)))


def main():
    root = o.HERE
    real = json.loads((root / "orion-results.json").read_text())
    ans = {}
    ans["orion"] = {}
    for key in ("single", "double", "triple"):
        for i, w in enumerate(real[key]["pair"]):
            state = np.array(w["state_eci_km_km_s"])
            # Fresh Python DOPRI5 OEL propagator; integrate to the emission epoch
            # and use finite Earth barycentric displacement directly for each ray.
            alt = OrbitPropagator(
                integrator="dopri5", plugins=o.engine.plugins, adaptive_atol=1e-11, adaptive_rtol=5e-13
            )
            tnodes = np.unique(np.r_[0, np.arange(0, o.times[-1], 300), o.times])
            history = []
            last = 0.0
            x = state.copy()
            for t in tnodes:
                if t > last:
                    x = alt.propagate(x, t - last, last, o.zero, o.env, o.ctx)
                last = t
                history.append(x)
            history = np.array(history)
            bodies = [o.hpop_de440_positions_km(o.JD + t / o.DAY, o.ephem) for t in tnodes]
            moonpos = np.array([b["moon"] for b in bodies])
            speed = np.linalg.norm(history[:, 3:], axis=1)
            earthdist = np.linalg.norm(history[:, :3], axis=1)
            moondist = np.linalg.norm(history[:, :3] - moonpos, axis=1)
            received = []
            for j, t in enumerate(o.times):
                k = np.searchsorted(tnodes, t)
                lag = np.linalg.norm(history[k, :3] - o.sites[j]) / o.LIGHT_SPEED
                # Advance from a preceding node, or from a backward subsecond
                # local Taylor step at the first reception only.
                for _ in range(2):
                    emit = t - lag
                    if emit < 0:
                        acc = alt.acceleration_at(
                            t_s=0, x_eci=state, command_accel_eci_km_s2=o.zero, env=o.env, ctx=o.ctx
                        )
                        pe = state[:3] + emit * state[3:] + 0.5 * emit**2 * acc
                    else:
                        prev = np.searchsorted(tnodes, emit) - 1
                        pe = alt.propagate(history[prev], emit - tnodes[prev], tnodes[prev], o.zero, o.env, o.ctx)[:3]
                    be = o.hpop_de440_positions_km(o.JD + emit / o.DAY, o.ephem)["earth"]
                    delta = pe - o.sites[j] + (be - bodies[k]["earth"])
                    lag = np.linalg.norm(delta) / o.LIGHT_SPEED
                received.append(o.unit(delta))
            received = np.array(received)
            residual = o.distance(received, o.obs)
            fitcount = real[key]["training_rows"]
            idx = np.r_[np.arange(fitcount), 20]
            mean_el = np.rad2deg(np.arcsin(np.sum(received[idx] * o.unit(o.sites[idx]), axis=1)))
            ans["orion"][key + str(i)] = dict(
                independent_fit_max_as=max(residual[:fitcount]),
                independent_first_holdout_los=received[20].tolist(),
                max_los_replay_difference_as=max(o.distance(received, o.sky(state, 25))),
                min_sampled_earth_radius_km=min(earthdist),
                max_sampled_earth_radius_km=max(earthdist),
                min_sampled_moon_radius_km=min(moondist),
                max_sampled_speed_km_s=max(speed),
                min_elevation_deg=min(mean_el),
                min_lunar_clearance_deg=min(
                    angle_visibility(received[idx], o.sites[idx], moonpos[np.searchsorted(tnodes, o.times[idx])])
                ),
            )
            assert max(residual[:fitcount]) < 2.0
    # Retrospective consistency check after all forecasts were fixed.
    seed = np.array(real["triple"]["baseline"]["state_eci_km_km_s"])
    opt = o.least_squares(
        lambda q: o.residual(o.sky(seed + o.correction(q), 25), 25).ravel(), np.zeros(6), jac="3-point", max_nfev=30
    )
    full = seed + o.correction(opt.x)
    err = o.distance(o.sky(full, 25), o.obs)
    ans["all_data_posthoc"] = dict(
        state=full.tolist(), max_as=max(err), rms_as=float(np.sqrt(np.mean(err**2))), coast_ok=o.admissible_coast(full)
    )
    (root / "checks.json").write_text(json.dumps(ans, indent=2) + "\n")


if __name__ == "__main__":
    main()
