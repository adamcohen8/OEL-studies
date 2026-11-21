"""Study flight software fed by a deterministic request and reply trace."""

from dataclasses import replace
from hashlib import sha256
from math import sin
from pathlib import Path

import numpy as np
from certificate import certify, guaranteed_remaining, stopping_distance
from evaluate import DRIFT, HANDOFF_S, MAX_REVERSE_S, approach

from sim.flight_software.contracts import (
    ActuatorCommand,
    DataValidity,
    DiagnosticTelemetry,
    DockingMechanismCommand,
    DockingPortMeasurement,
    FlightSoftwareOutput,
    FrameId,
    IdealOwnStateMeasurement,
    IdealWrenchCommand,
    MeasurementEvent,
    PacketId,
    StackIdentity,
    TelemetryField,
    ValidityInterval,
)
from sim.utils.quaternion import quaternion_to_dcm_bn

NS = 1_000_000_000
CHASER_PERIOD_NS = 31_250_000
RECEIVER_PERIOD_NS = 62_500_000
REQUEST = 2
TOPIC = "study.docking.contact_reservation.v1"


class LeaseStack:
    numeric_backend = "python"

    def __init__(self, role, policy, grant_ns, lease_ns, messages, lag_ns,
                 epoch_ns, rate_ppm, receiver_accel, measurement_bias, promise_break_ns=None,
                 navigation_initialization="ideal", tilt_bound_rad=0.0, resolve_wrench=False):
        if navigation_initialization != "ideal":
            raise ValueError("The test requires the OEL port measurement producer.")
        self.role = role
        self.policy = policy
        self.grant_ns = grant_ns
        self.lease_ns = lease_ns
        self.messages = [tuple(message) for message in messages]
        self.lag_ns = lag_ns
        self.epoch_ns = epoch_ns
        self.rate_ppm = rate_ppm
        self.receiver_accel = receiver_accel
        self.gap_bias, self.speed_bias = measurement_bias
        self.promise_break_ns = promise_break_ns
        self.tilt_bound_rad = tilt_bound_rad
        self.resolve_wrench = resolve_wrench
        digest = sha256()
        for name in ("lease_stacks.py", "certificate.py", "evaluate.py"):
            digest.update(Path(__file__).with_name(name).read_bytes())
        self.identity = StackIdentity(f"study.contact_reservation.{role}", "1.0", 1,
                                      digest.hexdigest(), False)

    def local_ns(self, time_ns):
        return self.epoch_ns + time_ns * (1_000_000 + self.rate_ppm) // 1_000_000

    def boot(self, event):
        self.boot_id = event.boot_id
        self.mode = "waiting"
        self.command_ns = None
        self.sent_local_ns = None
        self.receiver_grant_local_ns = None
        self.certificate = None
        self.eligible_remaining = None
        self.decision_gap = None
        self.decision_speed = None
        self.seen = set()
        self.ignored = 0

    def snapshot(self):
        raise ValueError("The study does not use checkpoints.")

    def restore(self, snapshot):
        raise ValueError("The study does not use checkpoints.")

    def shutdown(self, event):
        pass

    def receiver(self, time_ns):
        local = self.local_ns(time_ns)
        if time_ns >= self.grant_ns and self.receiver_grant_local_ns is None:
            self.receiver_grant_local_ns = local
        armed = (self.receiver_grant_local_ns is not None
                 and local < self.receiver_grant_local_ns + self.lease_ns
                 and (self.promise_break_ns is None or time_ns < self.promise_break_ns))
        self.mode = "armed" if armed else "unarmed"
        return armed, self.receiver_accel

    def read_grant(self, time_ns):
        grant = None
        for index, (arrival_ns, request, duration_ns, processing_ns) in enumerate(self.messages):
            if time_ns < arrival_ns or index in self.seen:
                continue
            self.seen.add(index)
            if request == REQUEST and grant is None and self.mode == "waiting":
                grant = duration_ns / NS, processing_ns / NS
            else:
                self.ignored += 1
        return grant

    def chaser(self, time_ns, port):
        if self.sent_local_ns is None:
            self.sent_local_ns = self.local_ns(time_ns)
        grant = self.read_grant(time_ns)
        if port.phase != "approach":
            self.mode = "contact"
            return False, 0.0
        if self.mode == "waiting":
            control_period = CHASER_PERIOD_NS / NS
            c = approach(port.offset_m[2] + self.gap_bias,
                         -port.velocity_m_s[2] + self.speed_bias)
            c = replace(c, gap_min=c.gap_min - 0.2 * sin(self.tilt_bound_rad))
            reserve = 0.002
            brake = 0.08
            safe_to_brake = c.gap_min >= reserve + stopping_distance(
                c.speed_max, c.lag_max, brake, c.disturbance)
            safe_to_wait = c.gap_min >= reserve + stopping_distance(
                c.speed_max, c.lag_max + control_period, brake, c.disturbance)
            if grant is not None and safe_to_brake:
                lease, processing = grant
                elapsed = (self.local_ns(time_ns) - self.sent_local_ns + 1) / NS
                if self.policy == "receipt_anchor":
                    remaining = lease / (1 + DRIFT)
                elif self.policy == "fixed_reverse":
                    remaining = lease / (1 + DRIFT) - MAX_REVERSE_S
                elif self.policy == "measured_round_trip":
                    remaining = guaranteed_remaining(lease, elapsed, DRIFT, DRIFT)
                else:
                    remaining = guaranteed_remaining(lease, elapsed, DRIFT, DRIFT, processing)
                self.eligible_remaining = remaining
                self.decision_gap = 0.5 * (c.gap_min + c.gap_max)
                self.decision_speed = 0.5 * (c.speed_min + c.speed_max)
                self.certificate = certify(c, remaining - HANDOFF_S)
                self.change_mode("commit" if self.certificate else "abort_window", time_ns)
            elif not safe_to_wait:
                self.change_mode("abort_margin", time_ns)

        active = self.command_ns is not None and time_ns >= self.command_ns + self.lag_ns
        if self.mode == "commit":
            return True, self.certificate.acceleration if active else 0.0
        return False, -0.08 if active and self.mode.startswith("abort") else 0.0

    def change_mode(self, mode, time_ns):
        self.mode = mode
        self.command_ns = time_ns

    def step(self, batch):
        now = batch.invocation_time
        time_ns = now.ticks * now.tick_period_ns
        if self.role == "receiver":
            armed, acceleration = self.receiver(time_ns)
            period_ns = RECEIVER_PERIOD_NS
        else:
            port = next((event.payload.payload for event in batch.events
                         if event.quality.validity is DataValidity.VALID
                         and isinstance(event.payload, MeasurementEvent)
                         and isinstance(event.payload.payload, DockingPortMeasurement)), None)
            if port is None:
                raise RuntimeError("Missing current port measurement.")
            armed, acceleration = self.chaser(time_ns, port)
            period_ns = CHASER_PERIOD_NS
        force = (0.0, 0.0, 100 * acceleration)
        if self.role == "chaser" and self.resolve_wrench:
            own = next(event.payload.payload for event in batch.events
                       if event.quality.validity is DataValidity.VALID
                       and isinstance(event.payload, MeasurementEvent)
                       and isinstance(event.payload.payload, IdealOwnStateMeasurement))
            attitude = quaternion_to_dcm_bn(own.attitude_quat_body_from_inertial)
            force = tuple(map(float, attitude @ (-100 * acceleration * np.asarray(port.receiver_axes_eci[2]))))
        expires = replace(now, ticks=now.ticks + 2 * period_ns // now.tick_period_ns)
        commands = tuple(
            ActuatorCommand(
                PacketId(f"{batch.satellite_id}/reservation", self.boot_id, 2 * batch.invocation_id + index),
                batch.satellite_id, device, now, ValidityInterval(now, expires),
                FrameId(f"OEL/ACTUATOR/{batch.satellite_id}/{device}", "frames-v1"), payload,
            )
            for index, (device, payload) in enumerate((
                ("wrench", IdealWrenchCommand(force, (0.0, 0.0, 0.0))),
                ("docking", DockingMechanismCommand(arm_capture=armed)),
            ))
        )
        values = {
            "mode": self.mode, "armed": armed, "clock_ns": self.local_ns(time_ns),
            "ignored_packets": self.ignored, "commanded_axial_acceleration_m_s2": acceleration,
        }
        optional = {
            "remaining_s": self.eligible_remaining,
            "decision_gap_m": self.decision_gap,
            "decision_speed_m_s": self.decision_speed,
        }
        if self.certificate:
            optional.update(earliest_s=self.certificate.earliest, latest_s=self.certificate.latest,
                            contact_speed_min_m_s=self.certificate.contact_speed_min,
                            contact_speed_max_m_s=self.certificate.contact_speed_max,
                            selected_acceleration_m_s2=self.certificate.acceleration)
        values.update({key: value for key, value in optional.items() if value is not None})
        return FlightSoftwareOutput(batch.satellite_id, batch.invocation_id, commands,
                                    (DiagnosticTelemetry(TOPIC, now, tuple(
                                        TelemetryField(key, value) for key, value in values.items()
                                    )),))
