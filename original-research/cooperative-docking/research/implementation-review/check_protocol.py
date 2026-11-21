"""Compare the study receiver and a separate grant-interval specification."""
from copy import copy
from itertools import product
from types import SimpleNamespace

from comparator_stack import ComparatorStack, Link, Request

ACTIONS = ("request", "renew", "repeat", "lost_renewal", "abort", "reset", "tick")


def check(template, actions):
    receiver = copy(template)
    receiver.link = Link()
    receiver.boot(SimpleNamespace(boot_id=1))
    original = {}
    stop = False
    for tick, action in enumerate(actions):
        now = 10 * tick
        if action in ("abort", "reset"):
            receiver.request_receiver_stop()
            stop = True
        packet = {"request": (1, 30), "renew": (2, 20), "lost_renewal": (2, 20),
                  "repeat": (1, 30)}.get(action)
        if action == "repeat" and 1 not in original:
            packet = None
        lost = action == "lost_renewal"
        receiver.profile = dict(template.profile, drop=(False, lost))
        expected_reply = None
        if packet:
            identity, duration = packet
            repeated = identity in original
            receiver.link.requests.append(Request(identity, now, duration))
            if not repeated and not stop:
                original[identity] = (now, duration)
            if identity in original and not lost:
                expected_reply = (identity, original[identity][1], repeated)
        armed, _ = receiver.receiver(now)
        assert armed == any(start <= now < start + duration for start, duration in original.values())
        assert {k: v.expiry_local_ns for k, v in receiver.ledger.items()} == {
            k: start + duration for k, (start, duration) in original.items()}
        replies = [(r.number, r.lease_ns, r.replayed) for r in receiver.link.replies]
        assert replies == ([] if expected_reply is None else [expected_reply])
        receiver.link.replies.clear()


def main():
    profile = {"turn_ns": 0, "caps_ns": (30, 20), "reverse_ns": (0, 0), "drop": (False, False)}
    template = ComparatorStack("receiver", "fixed_adaptive", "contract-check", profile,
                               0, 0, 0, 0., (0., 0.))
    for actions in product(ACTIONS, repeat=6):
        check(template, actions)
    # Lost initial reply plus a later duplicate must never grant a fresh fixed age.
    receiver = copy(template)
    receiver.link = Link()
    receiver.boot(SimpleNamespace(boot_id=1))
    receiver.profile = dict(profile, turn_ns=4, drop=(True, False))
    receiver.link.requests.append(Request(1, 0, 30))
    receiver.receiver(4)
    receiver.profile = profile
    receiver.link.requests.append(Request(1, 35, 30))
    assert not receiver.receiver(35)[0]
    reply = receiver.link.replies[0]
    assert (reply.lease_ns, reply.turnaround_ns, reply.replayed) == (30, 3, True)
    assert receiver.ledger[1].expiry_local_ns == 34
    chaser = copy(template)
    chaser.boot(SimpleNamespace(boot_id=1))
    chaser.link = receiver.link
    chaser.requests_sent[1] = 0
    chaser.receive(35)
    assert not chaser.grants_seen
    print(f"{len(ACTIONS)**6} implementation/specification traces passed; stale retransmission ignored.")


if __name__ == "__main__":
    main()
