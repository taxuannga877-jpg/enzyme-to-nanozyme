"""Helpers for formatting parsed PDB metadata for the Flask API."""

from nanozyme_mining.utils.constants import METAL_RESIDUE_NAMES

# Atom name to functional-group description.
_ATOM_FUNCTIONAL_GROUP = {
    "NE2": "Imidazole N\u03b5 (His)",
    "ND1": "Imidazole N\u03b4 (His)",
    "SG": "Thiol S (Cys)",
    "SD": "Thioether S (Met)",
    "OD1": "Carboxyl O (Asp)",
    "OD2": "Carboxyl O (Asp)",
    "OE1": "Carboxyl O (Glu)",
    "OE2": "Carboxyl O (Glu)",
    "OG": "Hydroxyl O (Ser)",
    "OG1": "Hydroxyl O (Thr)",
    "OH": "Phenol O (Tyr)",
    "O": "Backbone carbonyl O",
    "N": "Backbone amide N",
    "SE": "Selenol Se (Sec)",
}


def _format_metal_sites_from_extractor(metal_sites_extracted: list) -> list:
    """Convert PDBMetalExtractor sites into the frontend response shape."""
    if not metal_sites_extracted:
        return "N/A"

    merged_sites = {}
    for site in metal_sites_extracted:
        chain = site.metal_chain if hasattr(site, "metal_chain") else "N/A"
        res_num = site.metal_residue_id if hasattr(site, "metal_residue_id") else "N/A"
        metal_type = site.metal_type if hasattr(site, "metal_type") else "Unknown"
        key = (chain, res_num, metal_type)

        if key not in merged_sites:
            merged_sites[key] = {
                "metal_type": metal_type,
                "metal_name": site.metal_name if hasattr(site, "metal_name") else "Unknown",
                "chain": chain,
                "residue_number": res_num,
                "functional_role": getattr(site, "functional_role", "unknown"),
                "coordinating_residues": [],
                "seen_residues": set(),
            }

        if hasattr(site, "coordinating_residues"):
            for res in site.coordinating_residues:
                res_key = (res.get("chain"), res.get("residue_id"), res.get("atom_name"))
                if res_key not in merged_sites[key]["seen_residues"]:
                    merged_sites[key]["seen_residues"].add(res_key)
                    atom_name = res.get("atom_name", "-")
                    merged_sites[key]["coordinating_residues"].append(
                        {
                            "atom": atom_name,
                            "functional_group": _ATOM_FUNCTIONAL_GROUP.get(atom_name, atom_name),
                            "residue": res.get("residue_name", "-"),
                            "chain": res.get("chain", "-"),
                            "number": res.get("residue_id", "-"),
                            "distance": res.get("distance", "-"),
                        }
                    )

    formatted_sites = []
    for _key, site_data in merged_sites.items():
        coord_residues = site_data["coordinating_residues"]
        formatted_sites.append(
            {
                "metal_type": site_data["metal_type"],
                "metal_name": site_data["metal_name"],
                "chain": site_data["chain"],
                "residue_number": site_data["residue_number"],
                "functional_role": site_data["functional_role"],
                "coordinating_residues": coord_residues if coord_residues else "N/A",
                "coordination_number": len(coord_residues),
                "coordination_geometry": _infer_geometry(
                    len(coord_residues), site_data.get("metal_type", "")
                ),
            }
        )

    formatted_sites.sort(key=lambda x: (x["chain"], x["residue_number"]))
    return formatted_sites if formatted_sites else "N/A"


def _infer_geometry(coord_num: int, metal_type: str = "") -> str:
    """
    Infer coordination geometry from coordination number and metal type.

    Resolution order:
      1. Metal-specific override for the (metal, CN) pair when known
      2. CN-only default, preserving the previous behavior for unknown metals
    """
    metal = (metal_type or "").upper().strip()

    overrides = {
        ("CU", 4): "square_planar",
        ("NI", 4): "square_planar",
        ("PD", 4): "square_planar",
        ("PT", 4): "square_planar",
        ("AU", 4): "square_planar",
        ("AG", 2): "linear",
        ("CU", 2): "linear",
        ("FE", 5): "square_pyramidal",
    }
    if (metal, coord_num) in overrides:
        return overrides[(metal, coord_num)]

    geometries = {
        2: "linear",
        3: "trigonal_planar",
        4: "tetrahedral",
        5: "trigonal_bipyramidal",
        6: "octahedral",
        7: "pentagonal_bipyramidal",
        8: "square_antiprismatic",
    }
    return geometries.get(coord_num, "unknown")


def _extract_metal_sites(parsed: dict) -> list:
    """Extract metal-site records from parsed PDB metadata."""
    metal_sites = []
    links = parsed.get("links") or []
    hets = parsed.get("hets") or {}
    hetnams = parsed.get("hetnams") or {}

    for het_id, het_info in hets.items():
        if het_id.upper() in METAL_RESIDUE_NAMES:
            name = "N/A"
            if het_id in hetnams and hetnams[het_id]:
                name = hetnams[het_id][0].get("text", "N/A")

            coord_residues = []
            for link in links:
                if link.get("res1_name") == het_id or link.get("res2_name") == het_id:
                    coord_residues.append(
                        {
                            "atom": link.get("atom2_name") or link.get("atom1_name"),
                            "residue": link.get("res2_name") or link.get("res1_name"),
                            "chain": link.get("chain2") or link.get("chain1"),
                            "number": link.get("res2_num") or link.get("res1_num"),
                            "distance": link.get("length") or "N/A",
                        }
                    )

            metal_sites.append(
                {
                    "metal_type": het_id,
                    "metal_name": name,
                    "chain": het_info.get("chain") or "N/A",
                    "residue_number": het_info.get("seq_num") or "N/A",
                    "coordinating_residues": coord_residues if coord_residues else "N/A",
                }
            )

    return metal_sites if metal_sites else "N/A"


def _extract_active_sites(parsed: dict) -> list:
    """Extract SITE and REMARK 800 active-site records."""
    sites = parsed.get("sites") or []
    remark_800 = parsed.get("remark_800") or []

    site_descriptions = {}
    for remark in remark_800:
        if "site_id" in remark:
            site_descriptions[remark["site_id"]] = ""
        elif "description" in remark:
            for site_id in site_descriptions:
                if not site_descriptions[site_id]:
                    site_descriptions[site_id] = remark["description"]
                    break

    active_sites = []
    for site in sites:
        site_id = site.get("site_id", "")
        residues = site.get("residues") or []
        res_list = [
            f"{residue['residue_name']} {residue['chain']}{residue['residue_number']}"
            for residue in residues
        ]

        active_sites.append(
            {
                "site_id": site_id or "N/A",
                "residues": ", ".join(res_list) if res_list else "N/A",
                "description": site_descriptions.get(site_id, "N/A"),
            }
        )

    return active_sites if active_sites else "N/A"


def _extract_structure_info(parsed: dict) -> dict:
    """Extract secondary-structure, sequence, and mutation records."""
    helices = parsed.get("helices") or []
    sheets = parsed.get("sheets") or []
    helix_count = len(helices)
    sheet_count = len(sheets)

    ssbonds = parsed.get("ssbonds") or []
    ssbond_list = []
    for ssbond in ssbonds:
        ssbond_list.append(
            f"{ssbond.get('res1_name','')} {ssbond.get('chain1','')}{ssbond.get('res1_num','')} - "
            f"{ssbond.get('res2_name','')} {ssbond.get('chain2','')}{ssbond.get('res2_num','')} "
            f"({ssbond.get('length','N/A')} \u00c5)"
        )

    seqres = parsed.get("seqres") or {}
    seq_info = {}
    for chain, residues in seqres.items():
        seq_info[chain] = " ".join(residues) if residues else "N/A"

    seqadv = parsed.get("seqadv") or []
    mutations = []
    for adv in seqadv:
        if "MUTATION" in (adv.get("conflict") or "").upper():
            mutations.append(
                f"{adv.get('db_res','')} \u2192 {adv.get('res_name','')} "
                f"({adv.get('chain','')}{adv.get('seq_num','')})"
            )

    return {
        "secondary_structure": (
            f"{helix_count} helices, {sheet_count} sheets"
            if helix_count or sheet_count
            else "N/A"
        ),
        "disulfide_bonds": ssbond_list if ssbond_list else "N/A",
        "sequence": seq_info if seq_info else "N/A",
        "mutations": mutations if mutations else "N/A",
    }


def _build_pdb_info_response(
    parsed: dict,
    pdb_id: str,
    ec_number: str,
    pdb_file=None,
    metal_sites_extracted=None,
) -> dict:
    """Build the /api/get_pdb_full_info response payload."""
    header = parsed.get("header") or {}
    compnd = parsed.get("compnd") or {}
    source = parsed.get("source") or {}

    basic_info = {
        "pdb_id": header.get("pdb_id") or pdb_id,
        "title": parsed.get("title") or "N/A",
        "classification": header.get("classification") or "N/A",
        "date": header.get("date") or "N/A",
        "ec_number": compnd.get("EC") or ec_number or "N/A",
        "molecule": compnd.get("MOLECULE") or "N/A",
        "chains": compnd.get("CHAIN") or "N/A",
        "mutation": compnd.get("MUTATION") or "N/A",
        "resolution": parsed.get("resolution") or "N/A",
    }

    source_info = {
        "organism": source.get("ORGANISM_SCIENTIFIC") or "N/A",
        "organism_common": source.get("ORGANISM_COMMON") or "N/A",
        "gene": source.get("GENE") or "N/A",
        "expression_system": source.get("EXPRESSION_SYSTEM") or "N/A",
    }

    if metal_sites_extracted:
        metal_sites = _format_metal_sites_from_extractor(metal_sites_extracted)
    else:
        metal_sites = _extract_metal_sites(parsed)

    active_sites = _extract_active_sites(parsed)
    structure_info = _extract_structure_info(parsed)

    return {
        "basic_info": basic_info,
        "source_info": source_info,
        "metal_sites": metal_sites,
        "active_sites": active_sites,
        "structure_info": structure_info,
    }
