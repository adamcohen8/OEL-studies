# Prospective Mac validation protocol

Freeze the case inputs, source hashes, predicted outcomes and scoring tolerances
before starting any simulation or live timing trial. Draw a fresh seed for the
sealed batch. The controllers receive scenario inputs and incoming messages.
Expected outcomes and margins stay in the scoring manifest. Keep every case,
including unexpected outcomes and execution failures.

The timing cohort contains eight paired profiles. For each profile, vary only
receiver promise duration to put otherwise identical encounters 6–16 ms inside
and outside the predicted admission boundary. Vary gap, closing speed, receiver
motion, navigation bias, logical clock rates, reply time and actuation lag
across pairs. The outside decision is a refusal by this certificate. It does
not assert that every controller or disturbance realization makes such an
approach physically impossible.

The geometry cohort contains five matched profiles. For each, place two frozen
tilts 3–6 mrad below and above the angle at which a nominal leading ring point
would first beat the center-only earliest bound. Replay both angles with the
finite-radius correction and with the center-only comparator. Require corrected
certificates to contain every observed first-contact time. Test whether the
center-only comparator misses the higher-tilt contacts and contains the lower
ones. Check the measured ring lead against the separate r sin(tilt) prediction.

Count a predicted successful case as validated only if OEL reports soft capture,
first-contact time and speed lie inside that controller's certificate (except
the predeclared center-only misses), and both the frozen and realized timing
bounds retain 40 ms to physical expiry. Count a predicted refusal as validated
only when no contact occurs and sampled clearance stays at or above 2 mm.
Compare the event with a separate nominal scalar forecast to a predeclared
time tolerance of 1 ms. This forecast is a fixture prediction. All observed
physics comes from OEL. Compare the decision time to within 1 ns and geometry
lead to within 0.01 mm. A planned center-only miss is scored separately from
the corrected system. Stop and disclose any genuine unplanned mismatch.

Run twelve live two-process traces on the same Mac, in six pairs with reverse
delay on either side of the computed threshold. Include a close pair within
6 ms and interleave old messages on selected trials. Nominally place the trials
using the requested sleeps. Score actual decisions from measured timestamps,
since operating-system scheduling can move a trial across the threshold.
Independently check every accepted contact deadline plus handoff against the
receiver's recorded expiration. Also record any nominal-versus-measured changes.

The cases are prospective and outcome-blind within the local workflow. This is
not an independently administered blind trial. The study uses one Mac, emulated
clock rates, ideal-wrench models and a synthetic contact interface.
