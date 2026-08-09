#!/usr/bin/env python3
"""Dependency-light inspection of E2N's versioned reaction specifications."""

from __future__ import annotations

import json

from nanozyme_mining.design.physchem_knowledge import knowledge_version
from nanozyme_mining.design.substrate_catalog import list_reaction_tasks


def main() -> None:
    tasks = list_reaction_tasks()
    payload = {
        "physchem_knowledge_schema": knowledge_version(),
        "reaction_task_count": len(tasks),
        "reaction_tasks": [
            {
                "nanozyme_type": task.nanozyme_type,
                "task_id": task.task_id,
                "assay_context": task.assay,
                "calculation_method": task.calculation.barrier_method,
                "evidence_class": task.calculation.validation_level,
            }
            for task in tasks
        ],
        "claim_boundary": (
            "These entries define screening tasks and evidence classes; they do not "
            "constitute experimental validation or kinetic-rate predictions."
        ),
    }
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
