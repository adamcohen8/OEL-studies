# Historical SGP4 replay

Run `python replay.py` offline with `input.tle` and OEL 1.1.0. The delivered `states.csv` contains 13 native TEME states from the explicit Python backend, with elapsed time in minutes, position in km, and velocity in km/s. No coordinate transformation was applied.

The exact input passed checksum validation, the dispatcher selected `OGP-SGP4`, every propagation call succeeded, all state components were finite, the reverse-order run reproduced every value exactly, and the CSV passed a numeric round-trip check. Centered position differences gave a maximum velocity residual of 0.00122186 km/s at 4320 min. Halving the step from 2 s to 1 s changed the derivative by at most 4.63e-6 km/s. This checks kinematic consistency, not historical accuracy.

The compiled Rust code path, using the same parsed elements, differed by up to 7.89 m in position at 360 min and 6.26 mm/s in velocity at 3960 min. Both differences exceeded the chosen diagnostic thresholds of 0.01 m and 0.01 mm/s, respectively. The workspace supports a reproducible Python-backend trajectory, but backend equivalence at those tolerances is not supported. Agreement with the published vectors remains untested.

Prior familiarity: I recognized catalog object 00005 as a common historical SGP4 verification case and was familiar with SGP4 and TEME conventions. I did not use remembered or retrieved reference vectors.