#!/usr/bin/env python3
"""Measure OS-thread growth while an installed fido package runs."""

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time


POLL_SECONDS = 0.01
DEFAULT_TIMEOUT_SECONDS = 120.0


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("library", help="installed fido library")
    parser.add_argument("layer", type=str.upper, choices=("GNU", "INTEL"),
                        help="GNU or INTEL")
    parser.add_argument("output", type=Path, help="output JSON path")
    parser.add_argument("--stress", action="store_true",
                        help="disable MKL dynamic threading and allow nested OpenMP")
    parser.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT_SECONDS)
    parser.add_argument("--assert-bounded", action="store_true",
                        help="fail if GNU peak LIVE threads (excluding dead/zombie "
                             "states X/x and Z) exceed live baseline + cores - 1")
    args = parser.parse_args()

    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    if args.assert_bounded and args.layer != "GNU":
        parser.error("--assert-bounded applies only to GNU")
    return args


def task_snapshot(pid):
    task_directory = Path("/proc/{}/task".format(pid))
    task_ids = os.listdir(task_directory)
    tasks = {}
    vanished_task_ids = []
    for task_id in task_ids:
        try:
            stat = (task_directory / task_id / "stat").read_text(encoding="utf-8")
        except (FileNotFoundError, ProcessLookupError):
            vanished_task_ids.append(task_id)
            continue
        name_start = stat.index("(")
        name_end = stat.rindex(")")
        fields = stat[name_end + 2:].split()
        tasks[task_id] = {
            "name": stat[name_start + 1:name_end],
            "state": fields[0],
            "cpu_ticks": int(fields[11]) + int(fields[12]),
        }
    if str(pid) not in tasks or tasks[str(pid)]["state"] in ("X", "x", "Z"):
        raise RuntimeError("live main PID missing from task snapshot")
    return len(task_ids), tasks, vanished_task_ids


def live_task_count(tasks):
    return sum(details["state"] not in ("X", "x", "Z") for details in tasks.values())


def snapshot_metadata(tasks, elapsed, raw_count=None, vanished_task_ids=()):
    return {
        "elapsed_seconds": elapsed,
        "thread_count": len(tasks) if raw_count is None else raw_count,
        "live_thread_count": live_task_count(tasks),
        "vanished_task_ids": sorted(vanished_task_ids, key=int),
        "tasks": {task_id: tasks[task_id] for task_id in sorted(tasks, key=int)},
    }


def update_task_diagnostics(current, first_seen, last_seen, elapsed):
    for task_id, details in current.items():
        if task_id not in first_seen:
            first_seen[task_id] = {**details, "elapsed_seconds": elapsed}
        last_seen[task_id] = {**details, "elapsed_seconds": elapsed}


def mapped_libraries(pid):
    libraries = set()
    try:
        with open("/proc/{}/maps".format(pid), encoding="utf-8") as maps:
            for line in maps:
                fields = line.rstrip("\n").split(None, 5)
                if len(fields) != 6:
                    continue
                path = fields[5]
                name = Path(path).name
                if path.startswith("/") and re.search(r"\.so(?:\.[0-9.]+)?$", name):
                    libraries.add(path.removesuffix(" (deleted)"))
    except OSError:
        pass
    return libraries


def linkage_check(libraries, layer):
    names = [Path(path).name.lower() for path in libraries]
    has_gnu_thread = any("libmkl_gnu_thread" in name for name in names)
    has_intel_thread = any("libmkl_intel_thread" in name for name in names)
    has_gomp = any(re.search(r"libgomp(?:\.so|$)", name) for name in names)
    has_iomp = any("libiomp" in name for name in names)
    if layer == "GNU":
        selected = has_gnu_thread and has_gomp and not has_intel_thread and not has_iomp
        expected = "libmkl_gnu_thread plus libgomp, with no Intel threading runtime"
    else:
        selected = has_intel_thread and has_iomp and has_gomp and not has_gnu_thread
        expected = "libmkl_intel_thread plus libiomp and libgomp"
    opposite_absent = ((not has_intel_thread and not has_iomp)
                       if layer == "GNU" else not has_gnu_thread)
    return {
        "expected": expected,
        "has_libmkl_gnu_thread": has_gnu_thread,
        "has_libmkl_intel_thread": has_intel_thread,
        "has_libgomp": has_gomp,
        "has_libiomp": has_iomp,
        "opposite_threading_layer_absent": opposite_absent,
        "passed": selected,
    }


def wait_for_marker(marker, process, deadline):
    while not marker.exists():
        if process.poll() is not None:
            return False
        if time.monotonic() >= deadline:
            return False
        time.sleep(POLL_SECONDS)
    return True


def terminate_process(process):
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except (OSError, ProcessLookupError):
        try:
            process.terminate()
        except OSError:
            pass
    try:
        process.wait(timeout=5.0)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except (OSError, ProcessLookupError):
            try:
                process.kill()
            except OSError:
                pass
        process.wait(timeout=5.0)


def parse_result(path):
    values = {}
    if not path.exists():
        return values
    for line in path.read_text(encoding="utf-8").splitlines():
        key, separator, value = line.partition("=")
        if not separator:
            continue
        if value == "NA":
            values[key] = None
        elif value == "TRUE":
            values[key] = True
        elif value == "FALSE":
            values[key] = False
        else:
            try:
                values[key] = int(value)
            except ValueError:
                try:
                    values[key] = float(value)
                except ValueError:
                    values[key] = value
    return values


def nested_build_info(values, prefix):
    suffix = prefix + "."
    return {key[len(suffix):]: value for key, value in values.items()
            if key.startswith(suffix)}


def tail(text, limit=4000):
    if not text:
        return ""
    return text[-limit:]


def write_json(path, result):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")


def run(args):
    cores = len(os.sched_getaffinity(0))
    if cores < 1:
        raise RuntimeError("CPU affinity exposes no CPUs")
    rscript = os.environ.get("RSCRIPT") or shutil.which("Rscript")
    if not rscript:
        raise RuntimeError("Rscript was not found")
    workload = Path(__file__).with_name("thread-workload.R")

    env = os.environ.copy()
    env.update({
        "MKL_NUM_THREADS": str(cores),
        "OMP_NUM_THREADS": str(cores),
        "OPENBLAS_NUM_THREADS": "1",
        "BLIS_NUM_THREADS": "1",
    })
    env.pop("MKL_DYNAMIC", None)
    env.pop("OMP_MAX_ACTIVE_LEVELS", None)
    if args.stress:
        env["MKL_DYNAMIC"] = "FALSE"
        env["OMP_MAX_ACTIVE_LEVELS"] = "2"

    with tempfile.TemporaryDirectory(prefix="fido-thread-benchmark-") as temporary:
        temporary_path = Path(temporary)
        ready = temporary_path / "ready"
        start = temporary_path / "start"
        calls_started = temporary_path / "calls-started"
        finished = temporary_path / "finished"
        r_result = temporary_path / "result"
        command = [
            rscript, "--vanilla", str(workload), str(args.library),
            str(ready), str(start), str(calls_started), str(finished),
            str(r_result), str(cores),
        ]
        process = subprocess.Popen(
            command,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True,
        )
        deadline = time.monotonic() + args.timeout
        timed_out = False
        initial_threads = None
        peak_threads = None
        initial_live_threads = None
        peak_live_threads = None
        samples = 0
        measured_wall_seconds = None
        libraries = set()
        baseline_tasks = None
        baseline_vanished_task_ids = []
        first_seen_tasks = {}
        last_seen_tasks = {}
        new_task_ids = set()
        peak_task_snapshot = None
        peak_live_snapshot = None
        thread_sampling_error = None
        try:
            if not wait_for_marker(ready, process, deadline):
                timed_out = time.monotonic() >= deadline
            if not timed_out and process.poll() is None:
                initial_threads, baseline_tasks, baseline_vanished_task_ids = task_snapshot(process.pid)
                initial_live_threads = live_task_count(baseline_tasks)
                libraries.update(mapped_libraries(process.pid))
                baseline_time = time.monotonic()
                update_task_diagnostics(
                    baseline_tasks, first_seen_tasks, last_seen_tasks, 0.0
                )
                start.write_text("start\n", encoding="utf-8")
                if not wait_for_marker(calls_started, process, deadline):
                    timed_out = time.monotonic() >= deadline
            if not timed_out and calls_started.exists():
                monitor_start = time.monotonic()
                while not finished.exists() and process.poll() is None:
                    try:
                        current, current_tasks, vanished_task_ids = task_snapshot(process.pid)
                    except (OSError, ValueError, IndexError, RuntimeError) as error:
                        thread_sampling_error = str(error)
                        break
                    current_live = live_task_count(current_tasks)
                    libraries.update(mapped_libraries(process.pid))
                    new_peak = peak_threads is None or current > peak_threads
                    new_live_peak = peak_live_threads is None or current_live > peak_live_threads
                    peak_threads = current if peak_threads is None else max(peak_threads, current)
                    peak_live_threads = (current_live if peak_live_threads is None
                                         else max(peak_live_threads, current_live))
                    samples += 1
                    elapsed = time.monotonic() - baseline_time
                    update_task_diagnostics(
                        current_tasks, first_seen_tasks, last_seen_tasks, elapsed
                    )
                    new_task_ids.update(set(current_tasks) - set(baseline_tasks))
                    if new_peak:
                        peak_task_snapshot = snapshot_metadata(
                            current_tasks, elapsed, current, vanished_task_ids
                        )
                    if new_live_peak:
                        peak_live_snapshot = snapshot_metadata(
                            current_tasks, elapsed, current, vanished_task_ids
                        )
                    if time.monotonic() >= deadline:
                        timed_out = True
                        break
                    next_sample = monitor_start + samples * POLL_SECONDS
                    time.sleep(max(0.0, next_sample - time.monotonic()))
                measured_wall_seconds = time.monotonic() - monitor_start
                libraries.update(mapped_libraries(process.pid))
            if timed_out or thread_sampling_error is not None:
                terminate_process(process)
            else:
                remaining = max(0.0, deadline - time.monotonic())
                try:
                    process.wait(timeout=remaining)
                except subprocess.TimeoutExpired:
                    timed_out = True
                    terminate_process(process)
            stdout, stderr = process.communicate()
        except BaseException:
            terminate_process(process)
            process.communicate()
            raise

        values = parse_result(r_result)
        linkage = linkage_check(libraries, args.layer)
        process_exit = process.returncode
        thread_sampling_passed = (
            initial_threads is not None and peak_threads is not None and
            initial_live_threads is not None and peak_live_threads is not None and
            initial_live_threads > 0 and peak_live_threads > 0 and
            peak_live_snapshot is not None and
            thread_sampling_error is None and
            samples > 0 and measured_wall_seconds is not None
        )
        bounded_limit = (None if initial_live_threads is None
                         else initial_live_threads + cores - 1)
        bounded_passed = (peak_live_threads is not None and bounded_limit is not None
                          and peak_live_threads <= bounded_limit)
        assert_passed = True
        if args.assert_bounded:
            assert_passed = bounded_passed

        status = "ok"
        failures = []
        if timed_out:
            status = "failed"
            failures.append("timeout")
        if process_exit != 0:
            status = "failed"
            failures.append("Rscript exit status {}".format(process_exit))
        if not thread_sampling_passed:
            status = "failed"
            failures.append("required thread samples were not captured")
        if thread_sampling_error is not None:
            status = "failed"
            failures.append("thread snapshot failed: {}".format(thread_sampling_error))
        if not linkage["passed"]:
            status = "failed"
            failures.append("linked library verification failed")
        for check in ("finite_parallel", "finite_serial", "comparison_passed",
                      "mkl_restored", "actual_ncores_passed", "build_info_passed"):
            if values.get(check) is not True:
                status = "failed"
                failures.append("workload check failed: {}".format(check))
        if args.assert_bounded and not assert_passed:
            status = "failed"
            failures.append("bounded-thread assertion failed")
        if not values:
            status = "failed"
            failures.append("workload result was not written")

        result = {
            "schema_version": 2,
            "status": status,
            "failures": failures,
            "fido_library": str(args.library),
            "layer": args.layer,
            "mode": "stress" if args.stress else "default",
            "runner_cores": cores,
            "timeout_seconds": args.timeout,
            "process_exit": process_exit,
            "timed_out": timed_out,
            "initial_thread_count": initial_threads,
            "peak_thread_count": peak_threads,
            "initial_live_thread_count": initial_live_threads,
            "peak_live_thread_count": peak_live_threads,
            "allowed_peak_live_thread_count": bounded_limit,
            "thread_sample_count": samples,
            "thread_sampling_passed": thread_sampling_passed,
            "thread_sampling_error": thread_sampling_error,
            "thread_diagnostics": {
                "sampling_interval_seconds": POLL_SECONDS,
                "clock_ticks_per_second": os.sysconf("SC_CLK_TCK"),
                "baseline_tasks": baseline_tasks,
                "baseline_vanished_task_ids": baseline_vanished_task_ids,
                "raw_peak_snapshot": peak_task_snapshot,
                "peak_live_snapshot": peak_live_snapshot,
                "excluded_dead_states": ["X", "x", "Z"],
                "new_task_ids": sorted(new_task_ids, key=int),
                "per_task": {
                    task_id: {
                        "first": first_seen_tasks[task_id],
                        "last": last_seen_tasks[task_id],
                        "cpu_tick_delta": (
                            last_seen_tasks[task_id]["cpu_ticks"] -
                            first_seen_tasks[task_id]["cpu_ticks"]
                        ),
                    }
                    for task_id in sorted(first_seen_tasks, key=int)
                },
            },
            "measured_wall_seconds": measured_wall_seconds,
            "linked_libraries": sorted(libraries),
            "linked_library_verification": linkage,
            "assert_bounded_requested": args.assert_bounded,
            "assert_bounded_passed": assert_passed,
            "build_info_before": nested_build_info(values, "before"),
            "build_info_after_parallel": nested_build_info(values, "after_parallel"),
            "build_info_after": nested_build_info(values, "after"),
            "r_thread_environment": nested_build_info(values, "thread_env"),
            "workload": {
                key: value for key, value in values.items()
                if not key.startswith(("before.", "after_parallel.", "after.",
                                      "thread_env."))
            },
            "stdout_tail": tail(stdout),
            "stderr_tail": tail(stderr),
        }
        write_json(args.output, result)
        return 0 if status == "ok" else 1


def main():
    args = parse_args()
    try:
        return run(args)
    except BaseException as error:
        result = {
            "schema_version": 2,
            "status": "failed",
            "failures": [str(error)],
        }
        write_json(args.output, result)
        print(str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
