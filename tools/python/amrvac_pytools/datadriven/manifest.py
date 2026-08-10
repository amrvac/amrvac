"""Project-manifest helpers shared by DataConstrain and DataDriven workflows."""

from __future__ import print_function

import json
from pathlib import Path


SCHEMA_VERSION = 2


def read_project_manifest(path):
    """Read and non-destructively upgrade a project manifest in memory."""

    path = Path(path)
    if not path.exists():
        return {"schema_version": SCHEMA_VERSION, "workflows": {}}
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (TypeError, ValueError) as error:
        raise ValueError("invalid project manifest: {}".format(path)) from error
    if int(manifest.get("schema_version", 1)) >= SCHEMA_VERSION:
        manifest.setdefault("workflows", {})
        return manifest

    legacy = dict(manifest)
    workflow_name = str(legacy.get("workflow", "DataConstrain")).strip().lower()
    key = "data_constrain" if workflow_name == "dataconstrain" else workflow_name
    upgraded = {
        "schema_version": SCHEMA_VERSION,
        "input_dir": legacy.get("input_dir"),
        "project_dir": legacy.get("project_dir"),
        "paths": legacy.get("paths", {}),
        "region": legacy.get("region"),
        "relaxation_grid": legacy.get("relaxation_grid"),
        "evolution_grid": legacy.get("evolution_grid"),
        "workflows": {key: legacy},
    }
    return {key: value for key, value in upgraded.items() if value is not None}


def update_project_manifest(path, shared=None, workflow=None, workflow_updates=None):
    """Merge shared state and one workflow namespace without clobbering peers."""

    from .writers import ensure_output_dir, write_json

    path = Path(path)
    ensure_output_dir(path.parent)
    manifest = read_project_manifest(path)
    manifest["schema_version"] = SCHEMA_VERSION
    if shared:
        for key, value in shared.items():
            if value is not None:
                manifest[key] = value
    if workflow:
        workflows = manifest.setdefault("workflows", {})
        state = workflows.setdefault(str(workflow), {})
        if workflow_updates:
            state.update(workflow_updates)
    write_json(path, manifest)
    return manifest


__all__ = ["SCHEMA_VERSION", "read_project_manifest", "update_project_manifest"]
