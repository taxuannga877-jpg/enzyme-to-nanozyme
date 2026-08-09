"""
Comprehensive PDB Parser - Parse All PDB Record Types
=======================================================

Parses various PDB record types:
- HEADER: PDB ID, classification, date
- TITLE: Structure title
- COMPND: Compound information (EC number, molecule name, chains)
- SOURCE: Organism information
- REMARK 2: Resolution
- REMARK 620: Metal coordination angles
- REMARK 800: Active site descriptions
- DBREF: Database references (UniProt)
- SEQRES: Complete sequence
- SEQADV: Sequence differences/mutations
- FORMUL: Chemical formulas
- CRYST1: Crystallographic information
- SSBOND: Disulfide bonds
- LINK: Non-standard connections
- CONECT: Atom connectivity
- SITE: Active site annotations
- HET, HETNAM, HETSYN: Heteroatom information
- HELIX, SHEET: Secondary structure
- HETATM: Ligands and cofactors
"""

import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any


class PDBParser:
    """Comprehensive PDB file parser for all record types."""

    def __init__(self):
        """Initialize the PDB parser."""
        pass

    def parse_pdb_file(self, pdb_file: Path) -> Dict[str, Any]:
        """
        Parse a PDB file and extract all record types.

        Args:
            pdb_file: Path to PDB file

        Returns:
            Dictionary containing all parsed information
        """
        if not pdb_file.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")

        with open(pdb_file, "r", encoding="utf-8", errors="replace") as f:
            return self.parse_pdb_lines(f.readlines(), source_path=pdb_file)

    def parse_pdb_lines(
        self,
        lines: List[str],
        source_path: Optional[Path] = None,
    ) -> Dict[str, Any]:
        """
        Parse already-loaded PDB lines and extract all record types.

        source_path is retained for diagnostics and backwards-compatible
        callers, but this method does not reopen it.
        """
        result = {
            # Header information
            "header": None,
            "title": "",
            "compnd": {},
            "source": {},
            # Experimental information
            "resolution": None,
            "remark_620": [],  # Metal coordination angles
            "remark_800": [],  # Active site descriptions
            # Sequence information
            "dbref": [],
            "seqres": {},
            "seqadv": [],
            "formul": {},
            # Crystal information
            "cryst1": None,
            # Existing fields
            "ssbonds": [],
            "links": [],
            "conects": [],
            "sites": [],
            "hets": {},
            "hetnams": {},
            "hetsyns": {},
            "helices": [],
            "sheets": [],
            "ligands": []
        }

        # Multi-line record accumulators
        title_lines = []
        compnd_lines = []
        source_lines = []

        for line in lines:
            record_type = line[0:6].strip()

            # Header information
            if record_type == "HEADER":
                result["header"] = self._parse_header(line)

            elif record_type == "TITLE":
                title_lines.append(line)

            elif record_type == "COMPND":
                compnd_lines.append(line)

            elif record_type == "SOURCE":
                source_lines.append(line)

            # Experimental information
            elif record_type == "REMARK":
                remark_num = line[7:10].strip() if len(line) > 10 else ""
                if remark_num == "2":
                    resolution = self._parse_remark_2(line)
                    if resolution is not None:
                        result["resolution"] = resolution
                elif remark_num == "620":
                    remark_620 = self._parse_remark_620(line)
                    if remark_620:
                        result["remark_620"].append(remark_620)
                elif remark_num == "800":
                    remark_800 = self._parse_remark_800(line)
                    if remark_800:
                        result["remark_800"].append(remark_800)

            # Sequence information
            elif record_type == "DBREF":
                dbref = self._parse_dbref(line)
                if dbref:
                    result["dbref"].append(dbref)

            elif record_type == "SEQRES":
                seqres = self._parse_seqres(line)
                if seqres:
                    chain = seqres["chain"]
                    if chain not in result["seqres"]:
                        result["seqres"][chain] = []
                    result["seqres"][chain].extend(seqres["residues"])

            elif record_type == "SEQADV":
                seqadv = self._parse_seqadv(line)
                if seqadv:
                    result["seqadv"].append(seqadv)

            elif record_type == "FORMUL":
                formul = self._parse_formul(line)
                if formul:
                    het_id = formul["het_id"]
                    result["formul"][het_id] = formul

            # Crystal information
            elif record_type == "CRYST1":
                result["cryst1"] = self._parse_cryst1(line)

            # Existing record types
            elif record_type == "SSBOND":
                ssbond = self._parse_ssbond(line)
                if ssbond:
                    result["ssbonds"].append(ssbond)

            elif record_type == "LINK":
                link = self._parse_link(line)
                if link:
                    result["links"].append(link)

            elif record_type == "CONECT":
                conect = self._parse_conect(line)
                if conect:
                    result["conects"].append(conect)

            elif record_type == "SITE":
                site = self._parse_site(line)
                if site:
                    result["sites"].append(site)

            elif record_type == "HET":
                het = self._parse_het(line)
                if het:
                    het_id = het["het_id"]
                    result["hets"][het_id] = het

            elif record_type == "HETNAM":
                hetnam = self._parse_hetnam(line)
                if hetnam:
                    het_id = hetnam["het_id"]
                    if het_id not in result["hetnams"]:
                        result["hetnams"][het_id] = []
                    result["hetnams"][het_id].append(hetnam)

            elif record_type == "HETSYN":
                hetsyn = self._parse_hetsyn(line)
                if hetsyn:
                    het_id = hetsyn["het_id"]
                    if het_id not in result["hetsyns"]:
                        result["hetsyns"][het_id] = []
                    result["hetsyns"][het_id].append(hetsyn)

            elif record_type == "HELIX":
                helix = self._parse_helix(line)
                if helix:
                    result["helices"].append(helix)

            elif record_type == "SHEET":
                sheet = self._parse_sheet(line)
                if sheet:
                    result["sheets"].append(sheet)

            elif record_type == "HETATM":
                # Parse HETATM for ligands and cofactors
                ligand = self._parse_hetatm_ligand(line)
                if ligand:
                    result["ligands"].append(ligand)

        # Process multi-line records
        result["title"] = self._parse_title_lines(title_lines)
        result["compnd"] = self._parse_compnd_lines(compnd_lines)
        result["source"] = self._parse_source_lines(source_lines)

        # Combine HET information
        result["het_info"] = self._combine_het_info(
            result["hets"], result["hetnams"], result["hetsyns"]
        )

        # Fallback: If no SSBOND records found, try calculating from coordinates
        if not result["ssbonds"]:
            try:
                calculated_ssbonds = self._calculate_ssbonds_from_lines(lines)
                if calculated_ssbonds:
                    result["ssbonds"] = calculated_ssbonds
                    result["_calculated_ssbonds"] = True  # Flag to indicate calculation
            except Exception as e:
                # Silently fail if calculation doesn't work
                pass

        # Fallback: Identify non-standard residues from ATOM records
        if not result["ligands"]:
            try:
                calculated_ligands = self._identify_nonstandard_residues_from_lines(lines)
                if calculated_ligands:
                    result["ligands"].extend(calculated_ligands)
                    result["_calculated_ligands"] = True  # Flag to indicate calculation
            except Exception as e:
                # Silently fail if calculation doesn't work
                pass

        return result

    # ==================== New parsing methods ====================

    def _parse_header(self, line: str) -> Optional[Dict]:
        """Parse HEADER record."""
        try:
            return {
                "classification": line[10:50].strip(),
                "date": line[50:59].strip(),
                "pdb_id": line[62:66].strip()
            }
        except (ValueError, IndexError):
            return None

    def _parse_title_lines(self, lines: List[str]) -> str:
        """Parse multi-line TITLE records."""
        title_parts = []
        for line in lines:
            text = line[10:80].strip() if len(line) > 10 else ""
            if text:
                title_parts.append(text)
        return " ".join(title_parts)

    def _parse_compnd_lines(self, lines: List[str]) -> Dict:
        """Parse multi-line COMPND records."""
        compnd = {}
        current_key = None
        current_value = []

        for line in lines:
            text = line[10:80].strip() if len(line) > 10 else ""
            if not text:
                continue

            # Check if this line starts a new key
            if ":" in text:
                # Save previous key-value
                if current_key and current_value:
                    value = " ".join(current_value).rstrip(";").strip()
                    compnd[current_key] = value

                parts = text.split(":", 1)
                current_key = parts[0].strip()
                current_value = [parts[1].strip()] if len(parts) > 1 else []
            else:
                # Continuation of previous value
                current_value.append(text)

        # Save last key-value
        if current_key and current_value:
            value = " ".join(current_value).rstrip(";").strip()
            compnd[current_key] = value

        return compnd

    def _parse_source_lines(self, lines: List[str]) -> Dict:
        """Parse multi-line SOURCE records."""
        source = {}
        current_key = None
        current_value = []

        for line in lines:
            text = line[10:80].strip() if len(line) > 10 else ""
            if not text:
                continue

            if ":" in text:
                if current_key and current_value:
                    value = " ".join(current_value).rstrip(";").strip()
                    source[current_key] = value

                parts = text.split(":", 1)
                current_key = parts[0].strip()
                current_value = [parts[1].strip()] if len(parts) > 1 else []
            else:
                current_value.append(text)

        if current_key and current_value:
            value = " ".join(current_value).rstrip(";").strip()
            source[current_key] = value

        return source

    def _parse_remark_2(self, line: str) -> Optional[float]:
        """Parse REMARK 2 for resolution."""
        try:
            if "RESOLUTION" in line.upper():
                parts = line.split()
                for i, part in enumerate(parts):
                    if part.upper() == "RESOLUTION.":
                        if i + 1 < len(parts):
                            return float(parts[i + 1])
                    try:
                        val = float(part)
                        if 0.5 < val < 10.0:
                            return val
                    except ValueError:
                        continue
        except (ValueError, IndexError):
            pass
        return None

    def _parse_remark_620(self, line: str) -> Optional[Dict]:
        """Parse REMARK 620 for metal coordination."""
        try:
            text = line[11:].strip() if len(line) > 11 else ""
            if not text or text.startswith("COORDINATION"):
                return None
            return {"raw": text}
        except (ValueError, IndexError):
            return None

    def _parse_remark_800(self, line: str) -> Optional[Dict]:
        """Parse REMARK 800 for active site descriptions."""
        try:
            text = line[11:].strip() if len(line) > 11 else ""
            if not text:
                return None
            if "SITE_IDENTIFIER" in text:
                parts = text.split(":")
                if len(parts) >= 2:
                    return {"site_id": parts[1].strip()}
            elif "SITE_DESCRIPTION" in text:
                parts = text.split(":")
                if len(parts) >= 2:
                    return {"description": parts[1].strip()}
            return {"raw": text}
        except (ValueError, IndexError):
            return None

    def _parse_dbref(self, line: str) -> Optional[Dict]:
        """Parse DBREF record for database references."""
        try:
            if len(line) < 68:
                return None
            return {
                "pdb_id": line[7:11].strip(),
                "chain": line[12:13].strip(),
                "seq_begin": int(line[14:18].strip()) if line[14:18].strip() else None,
                "seq_end": int(line[20:24].strip()) if line[20:24].strip() else None,
                "database": line[26:32].strip(),
                "db_accession": line[33:41].strip(),
                "db_id": line[42:54].strip(),
                "db_seq_begin": int(line[55:60].strip()) if line[55:60].strip() else None,
                "db_seq_end": int(line[62:67].strip()) if line[62:67].strip() else None
            }
        except (ValueError, IndexError):
            return None

    def _parse_seqres(self, line: str) -> Optional[Dict]:
        """Parse SEQRES record for sequence."""
        try:
            if len(line) < 19:
                return None
            chain = line[11:12].strip()
            residues = line[19:70].split()
            return {"chain": chain, "residues": residues}
        except (ValueError, IndexError):
            return None

    def _parse_seqadv(self, line: str) -> Optional[Dict]:
        """Parse SEQADV record for mutations."""
        try:
            if len(line) < 49:
                return None
            return {
                "pdb_id": line[7:11].strip(),
                "res_name": line[12:15].strip(),
                "chain": line[16:17].strip(),
                "seq_num": int(line[18:22].strip()) if line[18:22].strip() else None,
                "database": line[24:28].strip(),
                "db_accession": line[29:38].strip(),
                "db_res": line[39:42].strip(),
                "db_seq": int(line[43:48].strip()) if line[43:48].strip() else None,
                "conflict": line[49:70].strip() if len(line) > 49 else ""
            }
        except (ValueError, IndexError):
            return None

    def _parse_formul(self, line: str) -> Optional[Dict]:
        """Parse FORMUL record for chemical formula."""
        try:
            if len(line) < 19:
                return None
            return {
                "comp_num": int(line[8:10].strip()) if line[8:10].strip() else None,
                "het_id": line[12:15].strip(),
                "formula": line[19:70].strip()
            }
        except (ValueError, IndexError):
            return None

    def _parse_cryst1(self, line: str) -> Optional[Dict]:
        """Parse CRYST1 record for crystal info."""
        try:
            if len(line) < 55:
                return None
            return {
                "a": float(line[6:15].strip()) if line[6:15].strip() else None,
                "b": float(line[15:24].strip()) if line[15:24].strip() else None,
                "c": float(line[24:33].strip()) if line[24:33].strip() else None,
                "alpha": float(line[33:40].strip()) if line[33:40].strip() else None,
                "beta": float(line[40:47].strip()) if line[40:47].strip() else None,
                "gamma": float(line[47:54].strip()) if line[47:54].strip() else None,
                "space_group": line[55:66].strip() if len(line) > 55 else "",
                "z": int(line[66:70].strip()) if len(line) > 66 and line[66:70].strip() else None
            }
        except (ValueError, IndexError):
            return None

    def _parse_ssbond(self, line: str) -> Optional[Dict]:
        """
        Parse SSBOND record (disulfide bond).

        Format:
        SSBOND   1 CYS A   23    CYS A  203                          1555   1555  2.03
        """
        try:
            if len(line) < 80:
                return None

            return {
                "serial": int(line[7:10].strip()) if line[7:10].strip() else None,
                "res1_name": line[11:14].strip(),
                "chain1": line[15:16].strip(),
                "res1_num": int(line[17:21].strip()) if line[17:21].strip() else None,
                "res2_name": line[25:28].strip(),
                "chain2": line[29:30].strip(),
                "res2_num": int(line[31:35].strip()) if line[31:35].strip() else None,
                "sym1": line[59:65].strip() if len(line) > 65 else "",
                "sym2": line[66:72].strip() if len(line) > 72 else "",
                "length": float(line[73:78].strip()) if len(line) > 78 and line[73:78].strip() else None
            }
        except (ValueError, IndexError):
            return None

    def _parse_link(self, line: str) -> Optional[Dict]:
        """
        Parse LINK record (non-standard connection).

        Format:
        LINK         SG   CYS A  23                 SG   CYS A 203     1555   1555  2.03
        """
        try:
            if len(line) < 80:
                return None

            return {
                "atom1_name": line[12:16].strip(),
                "res1_name": line[17:20].strip(),
                "chain1": line[21:22].strip(),
                "res1_num": int(line[22:26].strip()) if line[22:26].strip() else None,
                "atom2_name": line[42:46].strip(),
                "res2_name": line[47:50].strip(),
                "chain2": line[51:52].strip(),
                "res2_num": int(line[52:56].strip()) if line[52:56].strip() else None,
                "sym1": line[59:65].strip() if len(line) > 65 else "",
                "sym2": line[66:72].strip() if len(line) > 72 else "",
                "length": float(line[73:78].strip()) if len(line) > 78 and line[73:78].strip() else None
            }
        except (ValueError, IndexError):
            return None

    def _parse_conect(self, line: str) -> Optional[Dict]:
        """
        Parse CONECT record (atom connectivity).

        Format:
        CONECT 1179  746 1184 1195 1203
        """
        try:
            if len(line) < 11:
                return None

            atom_serial = int(line[6:11].strip())
            connected = []

            # CONECT can have up to 4 connected atoms
            for i in range(4):
                start = 11 + i * 5
                end = start + 5
                if len(line) > end:
                    atom_id = line[start:end].strip()
                    if atom_id:
                        try:
                            connected.append(int(atom_id))
                        except ValueError:
                            pass

            return {
                "atom_serial": atom_serial,
                "connected_atoms": connected
            }
        except (ValueError, IndexError):
            return None

    def _parse_site(self, line: str) -> Optional[Dict]:
        """
        Parse SITE record (active site annotation).

        Format:
        SITE    1 AC1 4 HIS A 146  HIS A  57  HIS A  87  HIS A 119
        """
        try:
            if len(line) < 22:
                return None

            site_id = line[11:14].strip()
            num_residues = int(line[15:17].strip()) if line[15:17].strip() else 0

            residues = []
            # SITE can have up to 4 residues per line
            for i in range(4):
                start = 18 + i * 11
                if len(line) > start + 10:
                    res_name = line[start:start+3].strip()
                    chain = line[start+4:start+5].strip()
                    res_num_str = line[start+5:start+10].strip()
                    if res_name and res_num_str:
                        try:
                            res_num = int(res_num_str)
                            residues.append({
                                "residue_name": res_name,
                                "chain": chain,
                                "residue_number": res_num
                            })
                        except ValueError:
                            pass

            return {
                "site_id": site_id,
                "num_residues": num_residues,
                "residues": residues
            }
        except (ValueError, IndexError):
            return None

    def _parse_het(self, line: str) -> Optional[Dict]:
        """
        Parse HET record (heteroatom information).

        Format:
        HET    HEM  A1547      10
        """
        try:
            if len(line) < 30:
                return None

            return {
                "het_id": line[7:10].strip(),
                "chain": line[12:13].strip(),
                "seq_num": int(line[13:17].strip()) if line[13:17].strip() else None,
                "num_atoms": int(line[20:25].strip()) if len(line) > 25 and line[20:25].strip() else None
            }
        except (ValueError, IndexError):
            return None

    def _parse_hetnam(self, line: str) -> Optional[Dict]:
        """
        Parse HETNAM record (heteroatom name).

        Format:
        HETNAM     HEM HEME
        """
        try:
            if len(line) < 13:
                return None

            het_id = line[11:14].strip()
            continuation = line[14:15].strip() if len(line) > 14 else ""
            text = line[15:70].strip() if len(line) > 15 else ""

            return {
                "het_id": het_id,
                "continuation": continuation,
                "text": text
            }
        except (ValueError, IndexError):
            return None

    def _parse_hetsyn(self, line: str) -> Optional[Dict]:
        """
        Parse HETSYN record (heteroatom synonym).

        Format:
        HETSYN     HEM HEME; PROTOHEME IX
        """
        try:
            if len(line) < 13:
                return None

            het_id = line[11:14].strip()
            continuation = line[14:15].strip() if len(line) > 14 else ""
            text = line[15:70].strip() if len(line) > 15 else ""

            return {
                "het_id": het_id,
                "continuation": continuation,
                "text": text
            }
        except (ValueError, IndexError):
            return None

    def _parse_helix(self, line: str) -> Optional[Dict]:
        """
        Parse HELIX record (secondary structure).

        Format:
        HELIX    1   1 PRO A    1  ALA A    5  1                                   5
        """
        try:
            if len(line) < 40:
                return None

            return {
                "serial": int(line[7:10].strip()) if line[7:10].strip() else None,
                "helix_id": line[11:14].strip(),
                "init_res_name": line[15:18].strip(),
                "init_chain": line[19:20].strip(),
                "init_seq_num": int(line[21:25].strip()) if line[21:25].strip() else None,
                "end_res_name": line[27:30].strip(),
                "end_chain": line[31:32].strip(),
                "end_seq_num": int(line[33:37].strip()) if line[33:37].strip() else None,
                "helix_class": int(line[38:40].strip()) if len(line) > 40 and line[38:40].strip() else None,
                "comment": line[40:70].strip() if len(line) > 40 else "",
                "length": int(line[71:76].strip()) if len(line) > 76 and line[71:76].strip() else None
            }
        except (ValueError, IndexError):
            return None

    def _parse_sheet(self, line: str) -> Optional[Dict]:
        """
        Parse SHEET record (secondary structure).

        Format:
        SHEET    1   A 4 PRO A  1  ALA A  5  0
        """
        try:
            if len(line) < 38:
                return None

            return {
                "strand": int(line[7:10].strip()) if line[7:10].strip() else None,
                "sheet_id": line[11:14].strip(),
                "num_strands": int(line[14:16].strip()) if line[14:16].strip() else None,
                "init_res_name": line[17:20].strip(),
                "init_chain": line[21:22].strip(),
                "init_seq_num": int(line[22:26].strip()) if line[22:26].strip() else None,
                "end_res_name": line[28:31].strip(),
                "end_chain": line[32:33].strip(),
                "end_seq_num": int(line[33:37].strip()) if line[33:37].strip() else None,
                "sense": int(line[38:40].strip()) if len(line) > 40 and line[38:40].strip() else None
            }
        except (ValueError, IndexError):
            return None

    def _parse_hetatm_ligand(self, line: str) -> Optional[Dict]:
        """
        Parse HETATM record for ligands and cofactors.

        Format:
        HETATM 1234  C1  HEM A1547      12.345  67.890  12.345  1.00 50.00           C
        """
        try:
            if len(line) < 54:
                return None

            # Check if it's a common ligand/cofactor (not a standard amino acid)
            res_name = line[17:20].strip()

            # 过滤掉水分子（HOH）
            if res_name.upper() == "HOH":
                return None

            # Common ligands and cofactors
            common_ligands = {
                "HEM", "HEA", "HEB", "HEC", "HEM",  # Heme groups
                "ZN", "FE", "CU", "MN", "MG", "CA",  # Metal ions
                "NAD", "NADP", "FAD", "FMN",  # Cofactors
                "ATP", "ADP", "AMP",  # Nucleotides
                "GOL", "SO4", "PO4",  # Common ligands
            }

            # Only include if it's a known ligand/cofactor or not a standard amino acid
            standard_aa = {
                "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
                "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
                "THR", "TRP", "TYR", "VAL", "SEC", "PYL"
            }

            if res_name in standard_aa:
                return None  # Skip standard amino acids

            atom_serial = int(line[6:11].strip())
            atom_name = line[12:16].strip()
            chain = line[21:22].strip()
            res_seq = int(line[22:26].strip()) if line[22:26].strip() else None
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            element = line[76:78].strip() if len(line) > 76 else ""

            return {
                "atom_serial": atom_serial,
                "atom_name": atom_name,
                "residue_name": res_name,
                "chain": chain,
                "residue_number": res_seq,
                "coordinates": [x, y, z],
                "element": element,
                "is_ligand": res_name not in standard_aa
            }
        except (ValueError, IndexError):
            return None

    def _combine_het_info(
        self,
        hets: Dict[str, Dict],
        hetnams: Dict[str, List[Dict]],
        hetsyns: Dict[str, List[Dict]]
    ) -> Dict[str, Dict]:
        """
        Combine HET, HETNAM, and HETSYN information.

        Args:
            hets: Dictionary of HET records
            hetnams: Dictionary of HETNAM records
            hetsyns: Dictionary of HETSYN records

        Returns:
            Combined heteroatom information
        """
        combined = {}

        for het_id, het_info in hets.items():
            combined[het_id] = {
                "het_id": het_id,
                "chain": het_info.get("chain", ""),
                "seq_num": het_info.get("seq_num"),
                "num_atoms": het_info.get("num_atoms"),
                "names": [],
                "synonyms": []
            }

            # Add names
            if het_id in hetnams:
                for hetnam in hetnams[het_id]:
                    text = hetnam.get("text", "")
                    if text:
                        combined[het_id]["names"].append(text)

            # Add synonyms
            if het_id in hetsyns:
                for hetsyn in hetsyns[het_id]:
                    text = hetsyn.get("text", "")
                    if text:
                        # Split by semicolon
                        synonyms = [s.strip() for s in text.split(";")]
                        combined[het_id]["synonyms"].extend(synonyms)

        return combined

    def extract_chemical_bonds(
        self,
        parsed_data: Dict[str, Any],
        pdb_file: Optional[Path] = None,
        metal_sites: Optional[List] = None
    ) -> Dict[str, List[Dict]]:
        """
        Extract and organize chemical bonds from parsed data.

        Args:
            parsed_data: Parsed PDB data
            pdb_file: Optional path to PDB file (for fallback calculations)
            metal_sites: Optional list of metal sites (for calculating metal coordination links)

        Returns:
            Dictionary with bond types as keys
        """
        ssbonds = parsed_data.get("ssbonds", [])
        links = parsed_data.get("links", [])
        conects = parsed_data.get("conects", [])

        # Fallback: If no links found and we have metal sites, try calculating metal coordination links
        if not links and pdb_file and metal_sites:
            try:
                calculated_links = self._calculate_links_from_coords(pdb_file, metal_sites)
                if calculated_links:
                    links.extend(calculated_links)
            except Exception as e:
                # Silently fail if calculation doesn't work
                pass

        return {
            "ssbonds": ssbonds,
            "links": links,
            "conects": conects
        }

    def extract_ligands_and_cofactors(
        self,
        parsed_data: Dict[str, Any],
        het_info: Optional[Dict[str, Dict]] = None
    ) -> List[Dict]:
        """
        Extract and organize ligands and cofactors.

        Args:
            parsed_data: Parsed PDB data
            het_info: Optional HET information dictionary

        Returns:
            List of ligand/cofactor dictionaries
        """
        ligands = parsed_data.get("ligands", [])

        # Group by residue name
        ligand_groups = {}
        for ligand in ligands:
            res_name = ligand.get("residue_name", "UNK")
            if res_name not in ligand_groups:
                ligand_groups[res_name] = {
                    "residue_name": res_name,
                    "chain": ligand.get("chain", ""),
                    "residue_number": ligand.get("residue_number"),
                    "atoms": [],
                    "het_info": None
                }

            ligand_groups[res_name]["atoms"].append({
                "atom_serial": ligand.get("atom_serial"),
                "atom_name": ligand.get("atom_name"),
                "element": ligand.get("element"),
                "coordinates": ligand.get("coordinates")
            })

        # Add HET information if available
        if het_info:
            for res_name, group in ligand_groups.items():
                if res_name in het_info:
                    group["het_info"] = het_info[res_name]

        return list(ligand_groups.values())

    def extract_secondary_structure(self, parsed_data: Dict[str, Any]) -> Dict[str, List[Dict]]:
        """
        Extract secondary structure information.

        Args:
            parsed_data: Parsed PDB data

        Returns:
            Dictionary with 'helix' and 'sheet' keys
        """
        return {
            "helix": parsed_data.get("helices", []),
            "sheet": parsed_data.get("sheets", [])
        }

    def _calculate_ssbonds_from_coords(self, pdb_file: Path) -> List[Dict]:
        """
        Calculate disulfide bonds from atomic coordinates.

        Finds CYS residues and calculates SG-SG distances.
        If distance < 2.5 Å, identifies as disulfide bond.

        Args:
            pdb_file: Path to PDB file

        Returns:
            List of disulfide bond dictionaries
        """
        with open(pdb_file, "r", encoding="utf-8", errors="replace") as f:
            return self._calculate_ssbonds_from_lines(f)

    def _calculate_ssbonds_from_lines(self, lines) -> List[Dict]:
        """Calculate disulfide bonds from already-loaded PDB lines."""
        ssbonds = []
        cys_sg_atoms = []  # List of (res_name, chain, res_id, x, y, z)

        # Parse all ATOM records to find CYS SG atoms
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    res_name = line[17:20].strip()
                    atom_name = line[12:16].strip()

                    # Only process CYS residues with SG atoms
                    if res_name == "CYS" and atom_name == "SG":
                        chain = line[21:22].strip()
                        res_id = int(line[22:26].strip()) if line[22:26].strip() else None
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())

                        cys_sg_atoms.append({
                            "res_name": res_name,
                            "chain": chain,
                            "res_id": res_id,
                            "coords": np.array([x, y, z])
                        })
                except (ValueError, IndexError):
                    continue

        # Calculate distances between all SG-SG pairs
        for i in range(len(cys_sg_atoms)):
            for j in range(i + 1, len(cys_sg_atoms)):
                atom1 = cys_sg_atoms[i]
                atom2 = cys_sg_atoms[j]

                # Calculate distance
                distance = float(np.linalg.norm(atom1["coords"] - atom2["coords"]))

                # Disulfide bond distance threshold: 2.5 Å
                # (Standard is ~2.0-2.1 Å, but we allow some tolerance for predicted structures)
                if distance < 2.5:
                    ssbond = {
                        "res1_name": atom1["res_name"],
                        "chain1": atom1["chain"],
                        "res1_num": atom1["res_id"],
                        "res2_name": atom2["res_name"],
                        "chain2": atom2["chain"],
                        "res2_num": atom2["res_id"],
                        "length": round(distance, 2),
                        "calculated_from_coords": True  # Flag to indicate this was calculated
                    }
                    ssbonds.append(ssbond)

        return ssbonds

    def _calculate_links_from_coords(
        self,
        pdb_file: Path,
        metal_sites: Optional[List] = None
    ) -> List[Dict]:
        """
        Calculate non-standard links from atomic coordinates.

        Currently focuses on metal coordination bonds.
        Can be extended for other non-standard connections.

        Args:
            pdb_file: Path to PDB file
            metal_sites: Optional list of metal sites (from metal extractor)

        Returns:
            List of link dictionaries
        """
        links = []

        # If metal sites are provided, extract metal coordination bonds
        if metal_sites:
            for metal_site in metal_sites:
                metal_coords = np.array(metal_site.get("metal_coords", []))
                if len(metal_coords) != 3:
                    continue

                coordinating_residues = metal_site.get("coordinating_residues", [])
                metal_type = metal_site.get("metal_type", "UNK")

                for coord_res in coordinating_residues:
                    coord_coords = coord_res.get("coordinates", [])
                    if len(coord_coords) != 3:
                        continue

                    distance = float(np.linalg.norm(metal_coords - np.array(coord_coords)))

                    # Metal coordination distance threshold: 3.5 Å
                    if distance < 3.5:
                        link = {
                            "atom1_name": metal_type,
                            "res1_name": metal_type,
                            "chain1": metal_site.get("metal_chain", ""),
                            "res1_num": metal_site.get("metal_residue_id"),
                            "atom2_name": coord_res.get("atom_name", ""),
                            "res2_name": coord_res.get("residue_name", ""),
                            "chain2": coord_res.get("chain", ""),
                            "res2_num": coord_res.get("residue_id"),
                            "length": round(distance, 2),
                            "calculated_from_coords": True
                        }
                        links.append(link)

        return links

    def _identify_nonstandard_residues(self, pdb_file: Path) -> List[Dict]:
        """
        Identify non-standard residues (ligands/cofactors) from ATOM records.

        Standard amino acids: ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY,
        HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL

        Args:
            pdb_file: Path to PDB file

        Returns:
            List of ligand/cofactor dictionaries
        """
        with open(pdb_file, "r", encoding="utf-8", errors="replace") as f:
            return self._identify_nonstandard_residues_from_lines(f)

    def _identify_nonstandard_residues_from_lines(self, lines) -> List[Dict]:
        """Identify non-standard residues from already-loaded PDB lines."""
        standard_amino_acids = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
            "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
            "THR", "TRP", "TYR", "VAL", "SEC", "PYL"  # Include selenocysteine and pyrrolysine
        }

        # Group atoms by residue
        residue_groups = {}  # (chain, res_id, res_name) -> list of atoms

        for line in lines:
            if line.startswith("ATOM"):
                try:
                    res_name = line[17:20].strip()
                    chain = line[21:22].strip()
                    res_id = int(line[22:26].strip()) if line[22:26].strip() else None

                    # Skip standard amino acids
                    if res_name in standard_amino_acids:
                        continue

                    # This is a non-standard residue (ligand/cofactor)
                    key = (chain, res_id, res_name)
                    if key not in residue_groups:
                        residue_groups[key] = []

                    atom_serial = int(line[6:11].strip())
                    atom_name = line[12:16].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    element = line[76:78].strip() if len(line) > 76 else ""

                    residue_groups[key].append({
                        "atom_serial": atom_serial,
                        "atom_name": atom_name,
                        "element": element,
                        "coordinates": [x, y, z]
                    })
                except (ValueError, IndexError):
                    continue

        # Convert to ligand format
        ligands = []
        for (chain, res_id, res_name), atoms in residue_groups.items():
            # Only include if it's a small molecule (atom count < 50)
            # This filters out large non-standard residues that might be errors
            if len(atoms) < 50:
                ligands.append({
                    "atom_serial": atoms[0]["atom_serial"],
                    "atom_name": atoms[0]["atom_name"],
                    "residue_name": res_name,
                    "chain": chain,
                    "residue_number": res_id,
                    "coordinates": atoms[0]["coordinates"],
                    "element": atoms[0]["element"],
                    "is_ligand": True,
                    "all_atoms": atoms,  # Include all atoms for this residue
                    "identified_from_atom": True  # Flag to indicate this was from ATOM record
                })

        return ligands
