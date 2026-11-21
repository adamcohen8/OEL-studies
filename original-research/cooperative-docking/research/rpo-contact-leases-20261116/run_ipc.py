"""Live timing replay in two OS processes using a local pipe and logical clocks."""

import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

from certificate import certify, guaranteed_remaining
from evaluate import DRIFT, HANDOFF_S, approach

NS = 1_000_000_000
REQUEST_ID = 2
REQUEST_STATE = approach(0.022, 0.040)


def local_ns(physical_ns, epoch_ns, rate_ppm):
    return epoch_ns + physical_ns * (1_000_000 + rate_ppm) // 1_000_000


def wait_until(physical_ns):
    delay = (physical_ns - time.monotonic_ns()) / NS
    if delay > 0:
        time.sleep(delay)


def receiver(conn, result, lease_s, reverse_s, mode):
    request = conn.recv()
    received_ns = time.monotonic_ns()
    time.sleep(0.10)
    issued_ns = time.monotonic_ns()
    rate_ppm, epoch = 1000, -240 * NS
    processing = local_ns(issued_ns, epoch, rate_ppm) - local_ns(received_ns, epoch, rate_ppm) - 1
    expiry_ns = issued_ns + int(lease_s * NS / (1 + rate_ppm / 1_000_000))
    send_ns = issued_ns + round(reverse_s * NS)
    if mode == "late":
        wait_until(expiry_ns)
        actual_disarm_ns = time.monotonic_ns()
    if mode != "drop":
        wait_until(send_ns)
        if mode == "reordered":
            conn.send((request - 1, 2 * NS, 0))
        conn.send((request, round(lease_s * NS), max(0, processing)))
    if mode != "late":
        wait_until(expiry_ns)
        actual_disarm_ns = time.monotonic_ns()
    result.put(("receiver", {"rx_ns": received_ns, "issued_ns": issued_ns,
                             "expiry_ns": expiry_ns, "disarm_ns": actual_disarm_ns}))
    conn.close()


def chaser(conn, result, wait_s):
    rate_ppm, epoch = -1000, 1112 * NS
    sent_ns = time.monotonic_ns()
    local_sent = local_ns(sent_ns, epoch, rate_ppm)
    time.sleep(0.02)
    conn.send(REQUEST_ID)
    deadline = sent_ns + round(wait_s * NS)
    ignored = 0
    reply = None
    while conn.poll(max(0.0, (deadline - time.monotonic_ns()) / NS)):
        try:
            request, duration_ns, processing_ns = conn.recv()
        except EOFError:
            break
        if request == REQUEST_ID:
            reply = duration_ns, processing_ns
            break
        ignored += 1
    decided_ns = time.monotonic_ns()
    info = {"tx_ns": sent_ns, "decision_ns": decided_ns, "ignored": ignored, "decision": "no_grant"}
    if reply is not None:
        elapsed = (local_ns(decided_ns, epoch, rate_ppm) - local_sent + 1) / NS
        duration, processing = (item / NS for item in reply)
        remaining = guaranteed_remaining(duration, elapsed, DRIFT, DRIFT, processing)
        cert = certify(REQUEST_STATE, remaining - HANDOFF_S)
        info.update(receiver_processing_s=processing, local_elapsed_s=elapsed,
                    certified_remaining_s=remaining, decision="commit" if cert else "reject")
        if cert:
            info["latest_contact_delay_s"] = cert.latest
    result.put(("chaser", info))
    conn.close()


def experiment(context, mode, lease, reverse, wait, expected):
    end1, end2 = context.Pipe()
    result = context.Queue()
    owner = context.Process(target=receiver, args=(end1, result, lease, reverse, mode))
    client = context.Process(target=chaser, args=(end2, result, wait))
    owner.start()
    client.start()
    end1.close()
    end2.close()
    try:
        rows = dict(result.get(timeout=5) for _ in range(2))
    finally:
        for participant in (owner, client):
            participant.join(timeout=2)
            if participant.is_alive():
                participant.terminate()
                participant.join()
        result.close()
    assert owner.exitcode == client.exitcode == 0
    c, r = rows["chaser"], rows["receiver"]
    assert c["decision"] == expected
    assert c["tx_ns"] <= r["rx_ns"] <= r["issued_ns"] < r["expiry_ns"] <= r["disarm_ns"]
    if expected == "commit":
        margin = (r["expiry_ns"] - c["decision_ns"]) / NS - c["latest_contact_delay_s"] - HANDOFF_S
        assert margin >= 0
        c["margin_to_actual_expiry_after_handoff_s"] = margin
    return {"trace": mode, **rows}


if __name__ == "__main__":
    ctx = mp.get_context("spawn")
    configurations = (("early", 0.82, 0.03, 1.5, "commit"),
                      ("reordered", 0.82, 0.045, 1.5, "commit"),
                      ("late", 0.43, 0.58, 1.5, "reject"),
                      ("drop", 0.42, 0.0, 0.45, "no_grant"))
    results = [experiment(ctx, *case) for case in configurations]
    path = Path(sys.argv[1])
    path.write_text(json.dumps(results, indent=2) + "\n")
    print([(row["trace"], row["chaser"]["decision"],
            row["chaser"].get("margin_to_actual_expiry_after_handoff_s")) for row in results])
