"""Paired OEL flight software with a causal, deterministic renewal link."""

from dataclasses import dataclass, replace
from hashlib import sha256
from math import sin
from pathlib import Path

from certificate import certify, guaranteed_remaining, stopping_distance
from evaluate import DRIFT, HANDOFF_S, MAX_REVERSE_S, approach
from lease_stacks import CHASER_PERIOD_NS, NS, LeaseStack

from sim.flight_software.contracts import DiagnosticTelemetry, StackIdentity, TelemetryField

TOPIC = "study.docking.adaptive_comparison.v1"
REQUEST_TRANSIT_NS = 31_250_000
FIRST_REPLY_TIMEOUT_NS = 625_000_000
WAIT_TIMEOUT_NS = 1_187_500_000
TARGET_SPEED_M_S = 0.032
BRAKE_M_S2 = 0.08
RESERVE_M = 0.002
RING_RADIUS_M = 0.2


@dataclass(frozen=True)
class Request:
    number: int
    arrived_ns: int
    lease_ns: int


@dataclass(frozen=True)
class Grant:
    number: int
    arrives_ns: int
    lease_ns: int
    turnaround_ns: int
    replayed: bool = False


@dataclass(frozen=True)
class Reservation:
    lease_ns: int
    turnaround_ns: int
    expiry_local_ns: int


class Link:
    """Transport mailboxes only. Peers cannot read each other's decisions."""

    def __init__(self):
        self.requests = []
        self.replies = []


_LINKS = {}


class ComparatorStack(LeaseStack):
    def __init__(self, role, policy, session, link_profile, lag_ns, epoch_ns, rate_ppm,
                 receiver_accel, measurement_bias, tilt_bound_rad=0.0, resolve_wrench=False,
                 navigation_initialization="ideal"):
        super().__init__(
            role, policy, grant_ns=0, lease_ns=0, messages=[], lag_ns=lag_ns,
            epoch_ns=epoch_ns, rate_ppm=rate_ppm, receiver_accel=receiver_accel,
            measurement_bias=measurement_bias, tilt_bound_rad=tilt_bound_rad,
            resolve_wrench=resolve_wrench, navigation_initialization=navigation_initialization,
        )
        content = sha256()
        for name in ("comparator_stack.py", "lease_stacks.py", "certificate.py", "evaluate.py"):
            content.update(Path(__file__).with_name(name).read_bytes())
        self.identity = StackIdentity(f"study.docking.adaptive_comparison.{role}", "1.1", 1,
                                      content.hexdigest(), False)
        self.profile = link_profile
        self.link = _LINKS.setdefault(session, Link())
        self.adaptive = policy != "fixed_one_shot"
        self.speed_management = policy != "fixed_renew_coast"

    def boot(self, event):
        super().boot(event)
        self.requests_sent = {}
        self.grants_seen = {}
        # Each stack owns one session. Keep original grants, including expired ones.
        if not hasattr(self, "ledger"):
            self.ledger = {}
            self.stop_requested = False
        self.issued = (0, 0)
        self.received = 0
        self.chosen = 0
        self.slowed = False
        self.first_declined = False
        self.pending = {}
        self.realized = 0.0

    def request_receiver_stop(self):
        """Decline new work and keep already-issued obligations until expiry."""
        self.stop_requested = True

    def receiver(self, time_ns):
        local = self.local_ns(time_ns)
        ready = [request for request in self.link.requests
                 if request.arrived_ns + self.profile["turn_ns"] <= time_ns]
        self.link.requests[:] = [request for request in self.link.requests if request not in ready]
        for request in ready:
            if request.number not in (1, 2):
                continue
            index = request.number - 1
            reservation = self.ledger.get(request.number)
            replayed = reservation is not None
            if not replayed:
                if self.stop_requested:
                    continue
                lease = min(request.lease_ns, self.profile["caps_ns"][index])
                turnaround = max(0, local - self.local_ns(request.arrived_ns) - 1)
                reservation = Reservation(lease, turnaround, local + lease)
                self.ledger[request.number] = reservation
                self.issued = (request.number, lease)
            if not self.profile["drop"][index]:
                self.link.replies.append(Grant(
                    request.number, time_ns + self.profile["reverse_ns"][index],
                    reservation.lease_ns, reservation.turnaround_ns, replayed,
                ))
        armed = any(local < grant.expiry_local_ns for grant in self.ledger.values())
        self.mode = "armed" if armed else "unarmed"
        return armed, self.receiver_accel

    def request(self, time_ns, number):
        self.requests_sent[number] = self.local_ns(time_ns)
        limit = self.profile["initial_request_ns"] if number == 1 else self.profile["renewal_request_ns"]
        self.link.requests.append(Request(number, time_ns + REQUEST_TRANSIT_NS, limit))

    def receive(self, time_ns):
        arrived = [reply for reply in self.link.replies if reply.arrives_ns <= time_ns]
        self.link.replies[:] = [reply for reply in self.link.replies if reply.arrives_ns > time_ns]
        for reply in arrived:
            # A retransmission has no original-issue-to-receipt delay guarantee.
            if reply.replayed:
                continue
            if reply.number in self.requests_sent and reply.number not in self.grants_seen:
                self.received = reply.number
                self.grants_seen[reply.number] = (reply, self.local_ns(time_ns))

    def remaining(self, grant, receipt_local, current_local):
        fixed = (grant.lease_ns / NS / (1 + DRIFT) - MAX_REVERSE_S
                 - (current_local - receipt_local + 1) / NS / (1 - DRIFT))
        if self.policy not in ("corrected_adaptive", "elapsed_adaptive"):
            return fixed
        corrected = guaranteed_remaining(
            grant.lease_ns / NS, (current_local - self.requests_sent[grant.number] + 1) / NS,
            DRIFT, DRIFT, grant.turnaround_ns / NS if self.policy == "corrected_adaptive" else 0.0,
        )
        return max(fixed, corrected)

    def options(self, state, time_ns):
        for number, (grant, receipt_local) in reversed(self.grants_seen.items()):
            remaining = self.remaining(grant, receipt_local, self.local_ns(time_ns))
            horizon = min(remaining - HANDOFF_S, 0.99 * state.speed_min / state.disturbance)
            solution = certify(state, horizon)
            if solution is not None:
                yield number, remaining, solution

    def command_after_lag(self, time_ns, acceleration):
        self.pending[time_ns + self.lag_ns] = acceleration

    def chaser(self, time_ns, port):
        # The software actuator samples on chaser ticks. Apply each overdue
        # setpoint in due-time order so the latest due value wins.
        for due in sorted(due for due in self.pending if due <= time_ns):
            self.realized = self.pending.pop(due)
        if port.phase != "approach":
            self.mode = "contact"
            return False, 0.0
        self.receive(time_ns)
        if self.mode == "commit":
            return True, self.realized
        if self.mode.startswith("abort"):
            return False, self.realized
        if not self.requests_sent:
            self.request(time_ns, 1)

        closing_speed = float(-port.velocity_m_s[2]) + self.speed_bias
        state = approach(float(port.offset_m[2]) + self.gap_bias, closing_speed)
        state = replace(state, gap_min=state.gap_min - RING_RADIUS_M * sin(self.tilt_bound_rad))
        coast = self.realized == 0 and all(command == 0 for command in self.pending.values())
        options = self.options(state, time_ns)
        safe_to_brake = state.gap_min >= RESERVE_M + stopping_distance(
            state.speed_max, state.lag_max, BRAKE_M_S2, state.disturbance)
        safe_to_wait = state.gap_min >= RESERVE_M + stopping_distance(
            state.speed_max, state.lag_max + CHASER_PERIOD_NS / NS, BRAKE_M_S2, state.disturbance)
        option = next(options, None) if coast and safe_to_brake else None
        if option is not None:
            self.chosen, self.eligible_remaining, self.certificate = option
            self.change_mode("commit", time_ns)
            self.command_after_lag(time_ns, self.certificate.acceleration)
            return True, self.realized

        initial_failed = 1 in self.grants_seen or time_ns >= FIRST_REPLY_TIMEOUT_NS
        if not safe_to_wait or time_ns >= WAIT_TIMEOUT_NS or (initial_failed and not self.adaptive):
            self.change_mode("abort_margin" if not safe_to_wait else "abort_window", time_ns)
            self.command_after_lag(time_ns, -BRAKE_M_S2)
            return False, self.realized

        if initial_failed and 2 not in self.requests_sent:
            self.first_declined = 1 in self.grants_seen
            self.request(time_ns, 2)
        order_accel = 0.0
        if 2 in self.requests_sent:
            self.mode = "renewing"
            if self.speed_management and 2 not in self.grants_seen:
                order_accel = max(-0.05, min(0.0, 1.8 * (TARGET_SPEED_M_S - closing_speed)))
                self.slowed |= order_accel < 0
        self.command_after_lag(time_ns, order_accel)
        return False, self.realized

    def step(self, batch):
        output = super().step(batch)
        metrics = {
            "requests": len(self.requests_sent),
            "issued_grant": self.issued[0], "issued_duration_ns": self.issued[1],
            "received_grant": self.received, "chosen_grant": self.chosen,
            "entered_slowdown": self.slowed, "first_declined_before_renewal": self.first_declined,
        }
        diagnostic = DiagnosticTelemetry(
            TOPIC, batch.invocation_time, tuple(TelemetryField(key, value) for key, value in metrics.items()),
        )
        return replace(output, telemetry=output.telemetry + (diagnostic,))
