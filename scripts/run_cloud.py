#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path

from pulser import QPUBackend, Sequence
from pulser.backend import BackendConfig
from pulser_pasqal import PasqalCloud


def required_environment() -> tuple[str, str, str]:
    names = ("PASQAL_USERNAME", "PASQAL_PASSWORD", "PASQAL_PROJECT_ID")
    missing = [name for name in names if not os.environ.get(name)]
    if missing:
        raise SystemExit(
            "Missing environment variables: " + ", ".join(missing)
        )
    return tuple(os.environ[name] for name in names)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate and optionally submit a sequence to Pasqal Cloud."
    )
    parser.add_argument("sequence", type=Path)
    parser.add_argument("--device", help="Cloud device name; defaults to first available QPU.")
    parser.add_argument("--shots", type=int, default=500)
    parser.add_argument(
        "--submit",
        action="store_true",
        help="Enable QPU submission. Without this flag no shots are spent.",
    )
    parser.add_argument(
        "--yes",
        action="store_true",
        help="Skip the final interactive QPU-spend confirmation.",
    )
    args = parser.parse_args()
    if not 1 <= args.shots <= 2000:
        raise SystemExit("--shots must be between 1 and 2000.")

    username, password, project_id = required_environment()
    connection = PasqalCloud(
        username=username,
        password=password,
        project_id=project_id,
    )
    devices = connection.fetch_available_devices()
    if not devices:
        raise SystemExit("Pasqal Cloud returned no available devices.")
    print("Available devices:", ", ".join(sorted(devices)))

    device_name = args.device or next(iter(devices))
    if device_name not in devices:
        raise SystemExit(f"Unknown device {device_name!r}.")
    device = devices[device_name]
    sequence = Sequence.from_abstract_repr(
        args.sequence.read_text(encoding="utf-8")
    ).with_new_device(device, strict=True)
    device.validate_register(sequence.register)
    print(
        f"Validated {args.sequence} for {device_name}: "
        f"{len(sequence.register.qubit_ids)} atoms, "
        f"{sequence.get_duration()} ns."
    )

    if not args.submit:
        print("Validation only. Re-run with --submit to request a QPU job.")
        return
    if not args.yes:
        answer = input(
            f"Spend {args.shots} hardware shots on {device_name}? Type YES: "
        )
        if answer != "YES":
            raise SystemExit("Submission cancelled.")

    backend = QPUBackend(
        sequence,
        connection,
        config=BackendConfig(default_num_shots=args.shots),
    )
    remote = backend.run(wait=False)
    print("Batch ID:", remote.batch_id)
    print("Job IDs:", ", ".join(map(str, remote.job_ids)))


if __name__ == "__main__":
    main()
