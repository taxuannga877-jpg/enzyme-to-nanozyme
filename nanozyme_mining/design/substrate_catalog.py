"""Activity-specific substrate and reaction-task catalog for catalysis screening."""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Dict, Iterable, List, Optional, Tuple

from .physchem_knowledge import get_activity_prototype, knowledge_version


@dataclass(frozen=True)
class SubstrateSpec:
    name: str
    smiles: str
    role: str
    copies: int = 1
    charge: int = 0
    spin: int = 1
    anchor_elements: Tuple[str, ...] = ("O", "N", "S", "P")
    target_distance: float = 2.8
    description: str = ""

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "smiles": self.smiles,
            "role": self.role,
            "copies": self.copies,
            "charge": self.charge,
            "spin": self.spin,
            "anchor_elements": list(self.anchor_elements),
            "target_distance": self.target_distance,
            "description": self.description,
        }


@dataclass(frozen=True)
class TransitionStateSpec:
    label: str
    kind: str
    coordinate: str
    substrate_names: Tuple[str, ...]
    reactive_bond_elements: Tuple[str, str] = ("", "")
    final_bond_distance: Optional[float] = None
    images: int = 7
    description: str = ""

    def to_dict(self) -> dict:
        return {
            "label": self.label,
            "kind": self.kind,
            "coordinate": self.coordinate,
            "substrate_names": list(self.substrate_names),
            "reactive_bond_elements": list(self.reactive_bond_elements),
            "final_bond_distance": self.final_bond_distance,
            "images": self.images,
            "description": self.description,
        }


@dataclass(frozen=True)
class CalculationProtocolSpec:
    """Scientific method requirements for one activity-screening task."""

    mechanism_family: str
    barrier_method: str
    requires_charge: bool = False
    requires_spin: bool = False
    neb_allowed: bool = False
    recommended_backends: Tuple[str, ...] = ("tblite", "pyscf")
    mace_heads: Tuple[str, ...] = ("omat_pbe", "omol", "oc20_usemppbe")
    validation_level: str = "screening_proxy"
    rationale: str = ""
    condition_id: Optional[str] = None
    ph_range: Tuple[float, float] = ()
    microstates: Tuple[str, ...] = ()
    solvent: str = "water"
    explicit_water_counts: Tuple[int, ...] = (0, 3, 6)
    coordinate_dimensions: int = 1
    knowledge_version: str = ""

    def to_dict(self) -> dict:
        return {
            "mechanism_family": self.mechanism_family,
            "barrier_method": self.barrier_method,
            "requires_charge": self.requires_charge,
            "requires_spin": self.requires_spin,
            "neb_allowed": self.neb_allowed,
            "recommended_backends": list(self.recommended_backends),
            "mace_heads": list(self.mace_heads),
            "validation_level": self.validation_level,
            "rationale": self.rationale,
            "condition_id": self.condition_id,
            "ph_range": list(self.ph_range),
            "microstates": list(self.microstates),
            "solvent": self.solvent,
            "explicit_water_counts": list(self.explicit_water_counts),
            "coordinate_dimensions": self.coordinate_dimensions,
            "knowledge_version": self.knowledge_version,
        }


@dataclass(frozen=True)
class ReactionTaskSpec:
    nanozyme_type: str
    task_id: str
    assay: str
    substrates: Tuple[SubstrateSpec, ...]
    transition_state: TransitionStateSpec
    calculation: CalculationProtocolSpec
    ml_task: str = "oc20"
    notes: Tuple[str, ...] = field(default_factory=tuple)

    def to_dict(self) -> dict:
        return {
            "nanozyme_type": self.nanozyme_type,
            "task_id": self.task_id,
            "assay": self.assay,
            "substrates": [s.to_dict() for s in self.substrates],
            "transition_state": self.transition_state.to_dict(),
            "calculation": self.calculation.to_dict(),
            "ml_task": self.ml_task,
            "notes": list(self.notes),
        }


_TMB = SubstrateSpec(
    name="TMB",
    smiles="Cc1cc(N)c(C)c(c1)c1c(C)c(N)cc(C)c1",
    role="chromogenic_electron_donor",
    anchor_elements=("N", "C"),
    target_distance=3.2,
    description="3,3',5,5'-tetramethylbenzidine assay donor.",
)

_H2O2 = SubstrateSpec(
    name="H2O2",
    smiles="OO",
    role="peroxide_oxidant",
    anchor_elements=("O",),
    target_distance=2.3,
)

_O2 = SubstrateSpec(
    name="O2",
    smiles="O=O",
    role="dioxygen_oxidant",
    spin=3,
    anchor_elements=("O",),
    target_distance=2.6,
)

_H2O = SubstrateSpec(
    name="H2O",
    smiles="O",
    role="water_nucleophile",
    anchor_elements=("O",),
    target_distance=2.4,
)


def _redox_protocol(mechanism_family: str, rationale: str) -> CalculationProtocolSpec:
    return CalculationProtocolSpec(
        mechanism_family=mechanism_family,
        barrier_method="electronic_state_scan",
        requires_charge=True,
        requires_spin=True,
        neb_allowed=False,
        rationale=rationale,
    )


def _hydrolysis_protocol(*, charged: bool = False, rationale: str) -> CalculationProtocolSpec:
    return CalculationProtocolSpec(
        mechanism_family="hydrolysis",
        barrier_method="coordinate_scan",
        requires_charge=charged,
        requires_spin=False,
        neb_allowed=False,
        rationale=rationale,
    )


REACTION_TASKS: Dict[str, ReactionTaskSpec] = {
    "Peroxidase": ReactionTaskSpec(
        nanozyme_type="Peroxidase",
        task_id="pod_tmb_h2o2",
        assay="TMB + H2O2 oxidation",
        substrates=(_TMB, _H2O2),
        transition_state=TransitionStateSpec(
            label="peroxide O-O activation",
            kind="bond_scission",
            coordinate="stretch H2O2 O-O while one O approaches the metal/oxo site",
            substrate_names=("H2O2",),
            reactive_bond_elements=("O", "O"),
            final_bond_distance=2.2,
            description="Initial NEB guess for peroxide activation before TMB oxidation.",
        ),
        calculation=_redox_protocol(
            "peroxide_redox",
            "Peroxide activation changes metal oxidation/spin state; geometry-only MACE NEB is not a kinetic barrier.",
        ),
        notes=(
            "Use adsorption energy for TMB and H2O2 separately and together.",
            "Barrier ranking should focus on H2O2 activation and TMB oxidation follow-up states.",
        ),
    ),
    "Catalase": ReactionTaskSpec(
        nanozyme_type="Catalase",
        task_id="cat_h2o2_disproportionation",
        assay="2 H2O2 disproportionation",
        substrates=(
            SubstrateSpec(
                name="H2O2",
                smiles="OO",
                role="peroxide_reductant_oxidant",
                copies=2,
                anchor_elements=("O",),
                target_distance=2.3,
            ),
        ),
        transition_state=TransitionStateSpec(
            label="peroxide O-O cleavage/proton transfer",
            kind="bond_scission",
            coordinate="activate one H2O2 O-O bond with a second H2O2/proton relay nearby",
            substrate_names=("H2O2",),
            reactive_bond_elements=("O", "O"),
            final_bond_distance=2.2,
        ),
        calculation=_redox_protocol(
            "peroxide_disproportionation",
            "Catalase turnover contains coupled peroxide redox and proton transfer steps.",
        ),
        notes=("Compare H2O2 binding that is strong enough to activate but not poison the site.",),
    ),
    "Oxidase": ReactionTaskSpec(
        nanozyme_type="Oxidase",
        task_id="oxd_tmb_o2",
        assay="TMB + O2 oxidation",
        substrates=(_TMB, _O2),
        transition_state=TransitionStateSpec(
            label="O2 activation",
            kind="electron_transfer",
            coordinate="bind O2 at the metal site and elongate O-O for superoxo/peroxo formation",
            substrate_names=("O2",),
            reactive_bond_elements=("O", "O"),
            final_bond_distance=1.6,
        ),
        calculation=_redox_protocol(
            "dioxygen_redox",
            "Triplet O2 activation and electron transfer require an explicit electronic state.",
        ),
    ),
    "Glucose Oxidase": ReactionTaskSpec(
        nanozyme_type="Glucose Oxidase",
        task_id="gox_glucose_o2",
        assay="glucose + O2 oxidation",
        substrates=(
            SubstrateSpec(
                name="D-glucose",
                smiles="OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O",
                role="sugar_substrate",
                anchor_elements=("O",),
                target_distance=2.8,
            ),
            _O2,
        ),
        transition_state=TransitionStateSpec(
            label="C-H/O2 hydride-transfer proxy",
            kind="hydride_transfer_proxy",
            coordinate="align glucose C1/OH region and O2 near the active metal for oxidative dehydrogenation",
            substrate_names=("D-glucose", "O2"),
            reactive_bond_elements=("C", "H"),
            final_bond_distance=1.6,
        ),
        calculation=_redox_protocol(
            "oxidative_dehydrogenation",
            "The coupled substrate oxidation/O2 reduction pathway cannot be represented by a single closed-shell bond stretch.",
        ),
    ),
    "Glutathione Peroxidase": ReactionTaskSpec(
        nanozyme_type="Glutathione Peroxidase",
        task_id="gpx_gsh_h2o2",
        assay="GSH + H2O2 reduction",
        substrates=(
            SubstrateSpec(
                name="GSH",
                smiles="N[C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(=O)O)C(=O)O",
                role="thiol_reductant",
                anchor_elements=("S", "O", "N"),
                target_distance=3.0,
            ),
            _H2O2,
        ),
        transition_state=TransitionStateSpec(
            label="thiol-peroxide proton/electron transfer",
            kind="proton_electron_transfer",
            coordinate="bring GSH sulfur and H2O2 oxygen into the metal second shell",
            substrate_names=("GSH", "H2O2"),
            reactive_bond_elements=("S", "O"),
            final_bond_distance=2.0,
        ),
        calculation=_redox_protocol(
            "proton_coupled_electron_transfer",
            "Selenol/thiol-peroxide chemistry requires proton and electron-state bookkeeping.",
        ),
    ),
    "Laccase": ReactionTaskSpec(
        nanozyme_type="Laccase",
        task_id="lac_catechol_o2",
        assay="catechol/ABTS-like phenol + O2 oxidation",
        substrates=(
            SubstrateSpec(
                name="catechol",
                smiles="Oc1ccccc1O",
                role="phenolic_electron_donor",
                anchor_elements=("O",),
                target_distance=3.0,
                description="Small ABTS-like phenolic model for fast screening.",
            ),
            _O2,
        ),
        transition_state=TransitionStateSpec(
            label="phenol electron transfer with O2 activation",
            kind="electron_transfer",
            coordinate="phenolic O near redox metal while O2 binds a second metal/edge site",
            substrate_names=("catechol", "O2"),
            reactive_bond_elements=("O", "O"),
            final_bond_distance=1.6,
        ),
        calculation=_redox_protocol(
            "multicopper_dioxygen_redox",
            "Laccase couples outer-sphere donor oxidation to multicopper O2 reduction with changing spin/oxidation states.",
        ),
    ),
    "Phosphatase": ReactionTaskSpec(
        nanozyme_type="Phosphatase",
        task_id="pho_pnpp_hydrolysis",
        assay="pNPP hydrolysis",
        substrates=(
            SubstrateSpec(
                name="pNPP",
                smiles="O=[N+]([O-])c1ccc(OP(=O)(O)O)cc1",
                role="phosphate_ester_substrate",
                anchor_elements=("P", "O"),
                target_distance=2.6,
                description="p-nitrophenyl phosphate.",
            ),
            _H2O,
        ),
        transition_state=TransitionStateSpec(
            label="P-O bond cleavage",
            kind="bond_scission",
            coordinate="stretch aryl P-O bond while water/nucleophile approaches phosphorus",
            substrate_names=("pNPP", "H2O"),
            reactive_bond_elements=("P", "O"),
            final_bond_distance=2.3,
            description="Fast NEB proxy for pNPP hydrolysis.",
        ),
        calculation=_hydrolysis_protocol(
            rationale="Use coupled P-O cleavage and nucleophilic O-P formation coordinates with explicit proton transfer validation.",
        ),
    ),
    "DNase": ReactionTaskSpec(
        nanozyme_type="DNase",
        task_id="dnase_dimethyl_phosphate_hydrolysis",
        assay="phosphodiester hydrolysis model",
        substrates=(
            SubstrateSpec(
                name="dimethyl phosphate",
                smiles="COP(=O)([O-])OC",
                role="phosphodiester_model",
                charge=-1,
                anchor_elements=("P", "O"),
                target_distance=2.6,
            ),
            _H2O,
        ),
        transition_state=TransitionStateSpec(
            label="phosphodiester P-O cleavage",
            kind="bond_scission",
            coordinate="stretch one P-O ester bond with water/nucleophile approaching P",
            substrate_names=("dimethyl phosphate", "H2O"),
            reactive_bond_elements=("P", "O"),
            final_bond_distance=2.4,
        ),
        calculation=_hydrolysis_protocol(
            charged=True,
            rationale="The phosphodiester model is anionic and requires charge-aware hydrolysis coordinates.",
        ),
    ),
    "Superoxide Dismutase": ReactionTaskSpec(
        nanozyme_type="Superoxide Dismutase",
        task_id="sod_superoxide_disproportionation",
        assay="superoxide disproportionation",
        substrates=(
            SubstrateSpec(
                name="superoxide",
                smiles="[O-][O]",
                role="superoxide_radical",
                copies=2,
                charge=-1,
                spin=2,
                anchor_elements=("O",),
                target_distance=2.4,
            ),
        ),
        transition_state=TransitionStateSpec(
            label="superoxide inner-sphere electron transfer",
            kind="electron_transfer",
            coordinate="bind superoxide at redox metal and compare first/second electron-transfer states",
            substrate_names=("superoxide",),
            reactive_bond_elements=("O", "O"),
            final_bond_distance=1.7,
        ),
        calculation=_redox_protocol(
            "redox_dismutation",
            "Two charged doublet superoxide states exchange electrons/protons through different metal oxidation states.",
        ),
    ),
    "Urease": ReactionTaskSpec(
        nanozyme_type="Urease",
        task_id="ure_urea_hydrolysis",
        assay="urea hydrolysis",
        substrates=(
            SubstrateSpec(
                name="urea",
                smiles="NC(N)=O",
                role="amide_substrate",
                anchor_elements=("O", "N"),
                target_distance=2.7,
            ),
            _H2O,
        ),
        transition_state=TransitionStateSpec(
            label="urea C-N cleavage",
            kind="nucleophilic_addition_elimination",
            coordinate="water/hydroxide attacks urea carbonyl carbon and weakens C-N bond",
            substrate_names=("urea", "H2O"),
            reactive_bond_elements=("C", "N"),
            final_bond_distance=2.0,
        ),
        calculation=_hydrolysis_protocol(
            rationale="Use coupled water/hydroxide attack and C-N cleavage coordinates; a pure C-N stretch is insufficient.",
        ),
    ),
}


ALIASES = {
    "POD": "Peroxidase",
    "CAT": "Catalase",
    "OXD": "Oxidase",
    "GOX": "Glucose Oxidase",
    "GPX": "Glutathione Peroxidase",
    "LAC": "Laccase",
    "PHO": "Phosphatase",
    "SOD": "Superoxide Dismutase",
    "URE": "Urease",
}


def normalize_activity(nanozyme_type: str) -> str:
    name = (nanozyme_type or "").strip()
    if name in REACTION_TASKS:
        return name
    upper = name.upper()
    if upper in ALIASES:
        return ALIASES[upper]
    for key in REACTION_TASKS:
        if key.lower() == name.lower():
            return key
    return name


def get_reaction_task(nanozyme_type: str) -> Optional[ReactionTaskSpec]:
    task = REACTION_TASKS.get(normalize_activity(nanozyme_type))
    if task is None:
        return None
    prototype = get_activity_prototype(task.nanozyme_type) or {}
    calculation = replace(
        task.calculation,
        condition_id=prototype.get("condition_id"),
        ph_range=tuple(float(value) for value in prototype.get("ph_range", ())),
        microstates=tuple(str(value) for value in prototype.get("microstates", ())),
        solvent="water",
        coordinate_dimensions=2 if task.calculation.mechanism_family == "hydrolysis" else 1,
        knowledge_version=knowledge_version(),
    )
    return replace(task, calculation=calculation)


def list_reaction_tasks(nanozyme_type: str = "") -> List[ReactionTaskSpec]:
    if nanozyme_type:
        task = get_reaction_task(nanozyme_type)
        return [task] if task else []
    return [get_reaction_task(k) for k in sorted(REACTION_TASKS)]


def expanded_substrates(task: ReactionTaskSpec) -> Iterable[SubstrateSpec]:
    for substrate in task.substrates:
        for _ in range(max(substrate.copies, 1)):
            yield substrate
