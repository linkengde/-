#!/usr/bin/env python3
"""Create portable-by-relink OVITO Pro scene states for corrected Round7 candidates."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import ovito
from ovito import scene
from ovito.io import import_file
from ovito.modifiers import (
    AssignColorModifier,
    ClearSelectionModifier,
    ComputePropertyModifier,
    ExpressionSelectionModifier,
    PythonScriptModifier,
    SelectTypeModifier,
)
from ovito.vis import Viewport

HERE = Path(__file__).resolve().parent
RUNS = HERE.parent
FULL_DUMP = RUNS / "round7_fiveframe_joint_analysis_to_1p6ps_20261109" / "round7_5frames_full_system_0p6_to_1p6ps.dump.gz"
JOINT_CSV = HERE / "round7_1p6ps_source_core_ag.csv"
CLUSTER_DUMP = HERE / "round7_1p6ps_largest_candidate_cluster.dump"


def ids_from_csv():
    with JOINT_CSV.open(newline="", encoding="utf-8") as stream:
        return [int(row["atom_id"]) for row in csv.DictReader(stream) if int(row["joint_candidate"]) == 1]


def set_view():
    viewport = scene.viewports.active_vp
    viewport.type = Viewport.Type.Perspective
    viewport.camera_dir = (1.0, -1.2, -0.8)
    viewport.camera_up = (0.0, 0.0, 1.0)
    viewport.zoom_all()


def save_full_system(candidate_ids):
    candidate_literal = repr(candidate_ids)
    script = f'''import numpy as np
def modify(frame, data):
    ids = np.asarray(data.particles["Particle Identifier"][:], dtype=np.int64)
    candidates = np.asarray({candidate_literal}, dtype=np.int64)
    values = np.isin(ids, candidates).astype(np.float64)
    if "joint_candidate" in data.particles.keys():
        prop = data.particles_.make_mutable(data.particles["joint_candidate"])
        prop[:] = values
    else:
        data.particles_.create_property("joint_candidate", data=values)
'''
    pipeline = import_file(str(FULL_DUMP))
    if pipeline.source.num_frames != 5:
        raise RuntimeError(f"Expected five frames, found {pipeline.source.num_frames}")
    pipeline.modifiers.append(PythonScriptModifier(script=script))
    pipeline.modifiers.append(SelectTypeModifier(types={1, 2, 3}))
    pipeline.modifiers.append(AssignColorModifier(color=(0.02, 0.02, 0.02)))
    pipeline.modifiers.append(SelectTypeModifier(types={4}))
    pipeline.modifiers.append(AssignColorModifier(color=(0.55, 0.75, 0.84)))
    pipeline.modifiers.append(ExpressionSelectionModifier(expression="joint_candidate == 1"))
    pipeline.modifiers.append(AssignColorModifier(color=(1.0, 0.40, 0.05)))
    pipeline.modifiers.append(ComputePropertyModifier(output_property="Transparency", expressions="0.38 * joint_candidate"))
    pipeline.modifiers.append(ClearSelectionModifier())
    pipeline.add_to_scene()
    scene.selected_pipeline = pipeline
    scene.anim.first_frame = 0
    scene.anim.last_frame = pipeline.source.num_frames - 1
    set_view()
    # Verify the saved scene's per-frame full-system and candidate counts.
    for frame in (0, 4):
        data = pipeline.compute(frame)
        count = int(np.count_nonzero(np.asarray(data.particles["joint_candidate"][:]) > 0.5))
        atom_count = len(data.particles.positions)
        if atom_count != 352637 or count != len(candidate_ids):
            raise RuntimeError(f"Frame {frame}: atoms={atom_count}, candidates={count}")
    scene.save(str(HERE / "round7_1p6ps_corrected_full_system.ovito"))


def save_largest_cluster():
    scene.pipelines.clear()
    pipeline = import_file(str(CLUSTER_DUMP))
    pipeline.modifiers.append(SelectTypeModifier(types={4}))
    pipeline.modifiers.append(AssignColorModifier(color=(1.0, 0.40, 0.05)))
    pipeline.modifiers.append(ClearSelectionModifier())
    pipeline.add_to_scene()
    scene.selected_pipeline = pipeline
    scene.anim.first_frame = 0
    scene.anim.last_frame = 0
    set_view()
    if len(pipeline.compute().particles.positions) != 693:
        raise RuntimeError("Largest-cluster OVITO pipeline does not contain 693 Ag atoms")
    scene.save(str(HERE / "round7_1p6ps_largest_cluster.ovito"))


if __name__ == "__main__":
    candidates = ids_from_csv()
    if len(candidates) != 874:
        raise RuntimeError(f"Expected 874 corrected joint candidates, found {len(candidates)}")
    save_full_system(candidates)
    save_largest_cluster()
    print("OVITO scenes saved and validated:", ovito.version_string)
