#!/usr/bin/env python3
"""Calculate CTA/CEA-2034 data curves from Klippel SPL Horizontal/Vertical exports.

Input files (tab-delimited Klippel exports, as e.g. provided by Audio Science Review speaker measurements):
    SPL Horizontal.txt
    SPL Vertical.txt

Output files are written to the requested CEA2034 and Early Reflections
subdirectories, plus Estimated In-Room Response.txt.

Dependency: numpy (pip install numpy)
No Spinorama repository or Spinorama modules are used.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np


# CTA-2034-A Appendix C, Table 7: spherical quadrangle weights at 10 degrees.
# The same value applies to +angle and -angle positions.
SOUND_POWER_WEIGHTS = {
    0: 0.000604486,
    10: 0.004730189,
    20: 0.008955027,
    30: 0.012387354,
    40: 0.014989611,
    50: 0.016868154,
    60: 0.018165962,
    70: 0.019006744,
    80: 0.019477787,
    90: 0.019629373,
    100: 0.019477787,
    110: 0.019006744,
    120: 0.018165962,
    130: 0.016868154,
    140: 0.014989611,
    150: 0.012387354,
    160: 0.008955027,
    170: 0.004730189,
    180: 0.000604486,
}


def parse_number(value: str) -> float:
    """Parse Klippel numeric cells, including thousands separators."""
    value = value.strip().replace("\u202f", "").replace("\xa0", "")
    if not value:
        return float("nan")
    # Klippel files use decimal points; commas in values such as 19,999.5
    # are thousands separators.
    return float(value.replace(",", ""))


def parse_angle(label: str) -> int:
    """Convert labels such as On-Axis, 10°, -20° into signed integer angles."""
    cleaned = label.strip().replace("−", "-")
    if cleaned.lower() in {"on-axis", "on axis", "0", "0°"}:
        return 0
    match = re.fullmatch(r"([+-]?\d+)\s*°?", cleaned)
    if not match:
        raise ValueError(f"Unrecognized angular column label: {label!r}")
    return int(match.group(1))


def load_klippel_file(filename: Path) -> Tuple[np.ndarray, Dict[int, np.ndarray]]:
    """Read a Klippel SPL Horizontal/Vertical tab-delimited export.

    Klippel's layout has two columns per angle: frequency and SPL. The first
    three rows are title/angle/column-description metadata.
    """
    with filename.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))

    if len(rows) < 4:
        raise ValueError(f"{filename.name}: file contains no measurement data")

    angle_row = rows[1]
    # The angle labels occupy columns 0, 2, 4, ...; intervening cells are blank.
    angle_columns: List[Tuple[int, int]] = []
    for column in range(0, len(angle_row), 2):
        label = angle_row[column].strip() if column < len(angle_row) else ""
        if not label:
            continue
        angle = parse_angle(label.strip('"'))
        angle_columns.append((angle, column))

    if not angle_columns:
        raise ValueError(f"{filename.name}: could not find angular column labels")

    frequency_rows: List[List[float]] = []
    for line_number, row in enumerate(rows[3:], start=4):
        if not row or not row[0].strip():
            continue
        try:
            values = [parse_number(cell) for cell in row]
        except ValueError as exc:
            raise ValueError(
                f"{filename.name}: invalid numeric value on line {line_number}: {exc}"
            ) from exc
        if len(values) < max(col for _, col in angle_columns) + 2:
            raise ValueError(
                f"{filename.name}: incomplete measurement row on line {line_number}"
            )
        frequency_rows.append(values)

    if not frequency_rows:
        raise ValueError(f"{filename.name}: no numeric measurement rows found")

    matrix = np.asarray(frequency_rows, dtype=float)
    curves: Dict[int, np.ndarray] = {}
    reference_frequency: np.ndarray | None = None

    for angle, column in angle_columns:
        this_frequency = matrix[:, column]
        spl = matrix[:, column + 1]
        if reference_frequency is None:
            reference_frequency = this_frequency
        elif not np.allclose(this_frequency, reference_frequency, rtol=0, atol=1e-6):
            raise ValueError(
                f"{filename.name}: frequency grid differs in the {angle} degree column"
            )
        if angle in curves:
            raise ValueError(f"{filename.name}: duplicate angle {angle} degrees")
        curves[angle] = spl

    assert reference_frequency is not None
    if np.any(~np.isfinite(reference_frequency)) or np.any(reference_frequency <= 0):
        raise ValueError(f"{filename.name}: invalid frequency values")
    if np.any(np.diff(reference_frequency) <= 0):
        raise ValueError(f"{filename.name}: frequencies are not strictly increasing")

    return reference_frequency, curves


def power_average(curves: Iterable[np.ndarray], weights: Iterable[float] | None = None) -> np.ndarray:
    """Energy/RMS average of SPL curves, returned in dB."""
    curve_list = [np.asarray(curve, dtype=float) for curve in curves]
    if not curve_list:
        raise ValueError("Cannot average an empty set of curves")

    values = np.vstack(curve_list)
    if np.any(~np.isfinite(values)):
        raise ValueError("A required measurement contains missing or non-finite SPL values")

    if weights is None:
        weight_array = np.ones(values.shape[0], dtype=float)
    else:
        weight_array = np.asarray(list(weights), dtype=float)
        if len(weight_array) != values.shape[0]:
            raise ValueError("Number of weights does not match number of curves")
        if np.any(weight_array < 0) or not np.any(weight_array > 0):
            raise ValueError("Weights must be non-negative and at least one must be positive")

    # SPL dB -> squared-pressure proportional quantity -> weighted average -> dB.
    power = np.power(10.0, values / 10.0)
    mean_power = np.average(power, axis=0, weights=weight_array)
    return 10.0 * np.log10(mean_power)


def require_angles(curves: Dict[int, np.ndarray], angles: Iterable[int], plane: str) -> List[np.ndarray]:
    """Return requested curves, with a useful error if an angle is missing."""
    result = []
    missing = []
    for angle in angles:
        if angle not in curves:
            missing.append(angle)
        else:
            result.append(curves[angle])
    if missing:
        raise ValueError(f"{plane} file is missing required angles: {missing}")
    return result


def write_curve(filename: Path, frequency: np.ndarray, spl: np.ndarray, title: str) -> None:
    """Write a two-column frequency/SPL file compatible with Spinorama text data."""
    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f'Freq "{title}"\n')
        for f, value in zip(frequency, spl):
            handle.write(f"{f:.6f}\t{value:.6f}\n")


def calculate(input_dir: Path, output_dir: Path) -> None:
    horizontal_path = input_dir / "SPL Horizontal.txt"
    vertical_path = input_dir / "SPL Vertical.txt"

    if not horizontal_path.is_file():
        raise FileNotFoundError(f"Input file not found: {horizontal_path}")
    if not vertical_path.is_file():
        raise FileNotFoundError(f"Input file not found: {vertical_path}")

    print("Reading horizontal measurements...")
    freq_h, horizontal = load_klippel_file(horizontal_path)
    print("Reading vertical measurements...")
    freq_v, vertical = load_klippel_file(vertical_path)

    if not np.allclose(freq_h, freq_v, rtol=0, atol=1e-6):
        raise ValueError("Horizontal and vertical files have different frequency grids")
    frequency = freq_h

    # The two measurement orbits share the 0-degree and 180-degree positions.
    # Use the horizontal orbit for these two shared points, as required to avoid
    # counting either shared spherical location twice in Sound Power.
    on_axis_h = horizontal[0]
    on_axis_v = vertical[0]
    if not np.allclose(on_axis_h, on_axis_v, rtol=0, atol=0.05, equal_nan=True):
        print("WARNING: horizontal and vertical on-axis curves differ by more than 0.05 dB; using horizontal.")

    # CTA-2034 Listening Window: on-axis, +/-10 degrees vertical, and +/-10/20/30 horizontal.
    listening_window = power_average(
        [on_axis_h]
        + require_angles(vertical, [-10, 10], "Vertical")
        + require_angles(horizontal, [-30, -20, -10, 10, 20, 30], "Horizontal")
    )

    # Individual early-reflection paths, as defined in CTA-2034-A section 5.2.
    floor_bounce = power_average(require_angles(vertical, [-20, -30, -40], "Vertical"))
    ceiling_bounce = power_average(require_angles(vertical, [40, 50, 60], "Vertical"))
    front_wall = power_average(require_angles(horizontal, [0, -10, 10, -20, 20, -30, 30], "Horizontal"))
    side_wall = power_average(require_angles(horizontal, [-40, 40, -50, 50, -60, 60, -70, 70, -80, 80], "Horizontal"))
    # Spinorama.org's corrected rear-wall calculation averages the complete
    # rear horizontal hemisphere: +/-90 through +/-170 degrees, plus 180.
    rear_wall_angles = [
        -170, -160, -150, -140, -130, -120, -110, -100, -90,
        90, 100, 110, 120, 130, 140, 150, 160, 170, 180,
    ]
    rear_wall = power_average(
        require_angles(horizontal, rear_wall_angles, "Horizontal")
    )

    # The total Early Reflections curve is the power average of all 26 directions
    # in the five CTA-2034 reflection groups. Each specified direction has equal
    # weight in this composite curve; the five separately exported components
    # above are also available for inspection.
    all_reflection_curves = (
        require_angles(vertical, [-20, -30, -40, 40, 50, 60], "Vertical")
        + require_angles(horizontal,
                         [0, -10, 10, -20, 20, -30, 30,
                          -40, 40, -50, 50, -60, 60, -70, 70, -80, 80,
                          -90, 90, 180], "Horizontal")
    )
    early_reflections = power_average(all_reflection_curves)

    # CTA-2034 Appendix C, Table 7 Sound Power weighting. There are 36 samples
    # in each orbit; the shared 0 and 180 degree positions are counted once.
    sound_power_curves: List[np.ndarray] = []
    sound_power_weights: List[float] = []
    for angle, curve in horizontal.items():
        sound_power_curves.append(curve)
        sound_power_weights.append(SOUND_POWER_WEIGHTS[abs(angle)])
    for angle, curve in vertical.items():
        if angle in (0, 180):
            continue  # already included from the horizontal orbit
        sound_power_curves.append(curve)
        sound_power_weights.append(SOUND_POWER_WEIGHTS[abs(angle)])

    if len(sound_power_curves) != 70:
        raise ValueError(
            f"Expected 70 unique polar measurements for Sound Power; found {len(sound_power_curves)}"
        )
    sound_power = power_average(sound_power_curves, sound_power_weights)

    # Directivity indices are defined as Listening Window minus the respective
    # curve. The calculated DI offset is added after it is determined below.

    # Estimated In-Room Response: 12% Listening Window, 44% Early Reflections,
    # 44% Sound Power. Weights apply to squared pressure, not dB values.
    estimated_in_room = power_average(
        [listening_window, early_reflections, sound_power],
        [0.12, 0.44, 0.44],
    )

    # User-requested sensitivity reference offset.
    sensitivity_mask = (frequency >= 300.0) & (frequency <= 3000.0)
    if not np.any(sensitivity_mask):
        raise ValueError("No frequency samples between 300 Hz and 3000 Hz")
    sensitivity = float(np.mean(on_axis_h[sensitivity_mask]))
    di_offset_value = math.floor((sensitivity - 40.0) / 10.0) * 10.0
    di_offset = np.full(frequency.shape, di_offset_value, dtype=float)

    # Apply the calculated offset to the DI curves, not to the underlying
    # Early Reflections or Sound Power response curves.
    early_reflections_di = listening_window - early_reflections + di_offset_value
    sound_power_di = listening_window - sound_power + di_offset_value

    cea_dir = output_dir / "CEA2034"
    reflections_dir = output_dir / "Early Reflections"

    outputs = [
        (cea_dir / "DI offset.txt", di_offset, "DI offset"),
        (cea_dir / "Early Reflections DI.txt", early_reflections_di, "Early Reflections DI"),
        (cea_dir / "Early Reflections.txt", early_reflections, "Early Reflections"),
        (cea_dir / "Listening Window.txt", listening_window, "Listening Window"),
        (cea_dir / "On Axis.txt", on_axis_h, "On Axis"),
        (cea_dir / "Sound Power DI.txt", sound_power_di, "Sound Power DI"),
        (cea_dir / "Sound Power.txt", sound_power, "Sound Power"),
        (reflections_dir / "Ceiling Bounce.txt", ceiling_bounce, "Ceiling Bounce"),
        (reflections_dir / "Floor Bounce.txt", floor_bounce, "Floor Bounce"),
        (reflections_dir / "Front Wall Bounce.txt", front_wall, "Front Wall Bounce"),
        (reflections_dir / "Rear Wall Bounce.txt", rear_wall, "Rear Wall Bounce"),
        (reflections_dir / "Side Wall Bounce.txt", side_wall, "Side Wall Bounce"),
        (reflections_dir / "Total Early Reflection.txt", early_reflections, "Total Early Reflection"),
        (output_dir / "Estimated In-Room Response.txt", estimated_in_room, "Estimated In-Room Response"),
    ]

    for filename, curve, title in outputs:
        write_curve(filename, frequency, curve, title)

    print("\nATTENTION:")
    print("This tool was generated by ChatGPT AI. To learn more about Spinorama data please refer to\nhttps://www.spinorama.org/ and https://github.com/spinorama-org/spinorama.\n")
    print("Please check the license terms under which the raw input data had been published. It may effect the calculated output\ndata as well.")
    license_path = input_dir / "LICENSE.txt"
    if license_path.is_file():
        print("I have found the following LICENSE.txt in the input directory:")
        with license_path.open() as f:
            print(f.read())
    print("\nCalculation complete.")
    print(f"Frequency points: {len(frequency)}")
    print(f"Frequency range: {frequency[0]:.4f} Hz to {frequency[-1]:.1f} Hz")
    print(f"Average on-axis sensitivity (300-3000 Hz): {sensitivity:.2f} dB")
    print(f"DI offset: {di_offset_value:.0f} dB")
    print(f"Output directory: {output_dir.resolve()}")
    print(f"Files written: {len(outputs)}")


def main() -> int:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Calculate CEA-2034 data curves from Klippel horizontal and vertical SPL exports.",
        epilog='''
Example:

  Download measurement ZIP archive from an ASR review (e.g. https://www.audiosciencereview.com/forum/index.php?threads/ascend-cbm-170se2-speaker-review.72683/),
  extract files and run

  CalcCEA2034.py 

  It creates

    CEA2034/Early Reflections DI.txt
    CEA2034/On Axis.txt
    CEA2034/Early Reflections.txt
    CEA2034/Sound Power.txt
    CEA2034/DI offset.txt
    CEA2034/Listening Window.txt
    CEA2034/Sound Power DI.txt
    Estimated In-Room Response.txt
    Early Reflections/Front Wall Bounce.txt
    Early Reflections/Total Early Reflection.txt
    Early Reflections/Floor Bounce.txt
    Early Reflections/Rear Wall Bounce.txt
    Early Reflections/Ceiling Bounce.txt
    Early Reflections/Side Wall Bounce.txt

  from the raw measurement data ("SPL Horizontal.txt", "SPL Vertical.txt").
'''
    )
    parser.add_argument(
        "input_dir", nargs="?", default=".",
        help="Folder containing \"SPL Horizontal.txt\" and \"SPL Vertical.txt\" (default: current folder)."
    )
    parser.add_argument(
        "--output-dir", default=None,
        help="Output folder (default: same as input folder)."
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else input_dir

    try:
        calculate(input_dir, output_dir)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
