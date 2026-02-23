from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple


@dataclass(frozen=True)
class AtomRow:
    nr: int
    atype: str
    resnr: int
    residue: str
    atom: str
    cgnr: int
    charge: float
    mass: Optional[float] = None


@dataclass(frozen=True)
class BondRow:
    ai: int
    aj: int
    funct: int
    params: Tuple[str, ...] = ()


@dataclass(frozen=True)
class AngleRow:
    ai: int
    aj: int
    ak: int
    funct: int
    params: Tuple[str, ...] = ()


@dataclass(frozen=True)
class DihedralRow:
    ai: int
    aj: int
    ak: int
    al: int
    funct: int
    params: Tuple[str, ...] = ()

# ----------------------------
# Parsers per section
# ----------------------------

def parse_atoms(tokens: List[str]) -> AtomRow:
    # nr type resnr residue atom cgnr charge [mass]
    if len(tokens) < 7:
        raise ValueError(f"[ atoms ] expects >=7 columns, got: {tokens}")
    mass = float(tokens[7]) if len(tokens) >= 8 else None
    return AtomRow(
        nr=int(tokens[0]),
        atype=tokens[1],
        resnr=int(tokens[2]),
        residue=tokens[3],
        atom=tokens[4],
        cgnr=int(tokens[5]),
        charge=float(tokens[6]),
        mass=mass,
    )

def int_or_str(x: str):
    try:
        return int(x)
    except ValueError:
        return x

def parse_bonds(tokens: List[str]) -> BondRow:
    if len(tokens) < 3:
        raise ValueError(f"[ bonds ] expects >=3 columns, got: {tokens}")
    return BondRow(
        ai=int(tokens[0]),
        aj=int(tokens[1]),
        funct=int_or_str(tokens[2]),
        params=tuple(tokens[3:]),
    )


def parse_angles(tokens: List[str]) -> AngleRow:
    if len(tokens) < 4:
        raise ValueError(f"[ angles ] expects >=4 columns, got: {tokens}")
    return AngleRow(
        ai=int(tokens[0]),
        aj=int(tokens[1]),
        ak=int(tokens[2]),
        funct=int_or_str(tokens[3]),
        params=tuple(tokens[4:]),
    )


def parse_dihedrals(tokens: List[str]) -> DihedralRow:
    if len(tokens) < 5:
        raise ValueError(f"[ dihedrals ] expects >=5 columns, got: {tokens}")
    return DihedralRow(
        ai=int(tokens[0]),
        aj=int(tokens[1]),
        ak=int(tokens[2]),
        al=int(tokens[3]),
        funct=int_or_str(tokens[4]),
        params=tuple(tokens[5:]),
    )


SECTION_PARSERS: Dict[str, Callable[[List[str]], Any]] = {
    "atoms": parse_atoms,
    "bonds": parse_bonds,
    "angles": parse_angles,
    "dihedrals": parse_dihedrals,
}


# ----------------------------
# MoleculeType object
# ----------------------------

@dataclass
class MoleculeType:
    name: str
    nrexcl: int
    sections: Dict[str, List[Any]] = field(default_factory=dict)

    def add_row(self, section: str, row: Any) -> None:
        self.sections.setdefault(section, []).append(row)

    def get_section(self, section: str) -> List[Any]:
        return self.sections.get(section.lower(), [])

    # typed convenience accessors
    @property
    def atoms(self) -> List[AtomRow]:
        return self.sections.get("atoms", [])  # type: ignore[return-value]

    @property
    def bonds(self) -> List[BondRow]:
        return self.sections.get("bonds", [])  # type: ignore[return-value]

    @property
    def angles(self) -> List[AngleRow]:
        return self.sections.get("angles", [])  # type: ignore[return-value]

    @property
    def dihedrals(self) -> List[DihedralRow]:
        return self.sections.get("dihedrals", [])  # type: ignore[return-value]


# ----------------------------

class ItpParseError(RuntimeError):
    pass


def strip_comment(line: str) -> str:
    return line.split(";", 1)[0].strip()


def is_section_header(line: str) -> bool:
    line = strip_comment(line)
    return line.startswith("[") and "]" in line


def parse_section_name(line: str) -> str:
    line = strip_comment(line)
    left = line.find("[")
    right = line.find("]", left + 1)
    if left == -1 or right == -1:
        raise ValueError(f"Invalid section header: {line!r}")
    return line[left + 1 : right].strip().lower()


def iter_lines(path: Path) -> Iterator[Tuple[int, str]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    for i, raw in enumerate(text.splitlines(), start=1):
        s = raw.strip()
        if not s:
            continue
        yield i, s


# ----------------------------
# Main: read .itp into MoleculeType objects
# ----------------------------

def read_itp_molecules(path: Path) -> List[MoleculeType]:
    molecules: List[MoleculeType] = []

    current_section: Optional[str] = None
    current_mol: Optional[MoleculeType] = None

    for lineno, raw in iter_lines(path):
        # record/ignore preprocessor directives for now
        if raw.startswith("#"):
            continue

        if is_section_header(raw):
            current_section = parse_section_name(raw)
            continue

        line = strip_comment(raw)
        if not line:
            continue

        if current_section is None:
            # fail fast: data before any section header
            raise ItpParseError(f"{path}:{lineno}: Data line outside any section: {raw!r}")

        tokens = line.split()

        if current_section == "moleculetype":
            if len(tokens) < 2:
                raise ItpParseError(f"{path}:{lineno}: [ moleculetype ] expects: name nrexcl")
            current_mol = MoleculeType(name=tokens[0], nrexcl=int(tokens[1]))
            molecules.append(current_mol)
            continue

        # If we see content for molecule sections before first moleculetype, decide policy:
        if current_mol is None:
            # Option 1: ignore global content
            continue

            # Option 2 (default here): error, because you said "Molecule contains all information for a single one"
            #raise ItpParseError(
            #    f"{path}:{lineno}: Section [{current_section}] appears before any [ moleculetype ]. "
            #    f"If your itp contains global parameter sections, we can support that explicitly."
            #)

        parser = SECTION_PARSERS.get(current_section)
        row = parser(tokens) if parser else list(tokens)  # unknown sections => keep raw tokens
        current_mol.add_row(current_section, row)

    return molecules

if __name__ == "__main__":
    mols = read_itp_molecules(Path("/home/guardian/Documents/TS2CG-v2.0/Tutorials/tut1/files/martini3/martini_v3.0_phospholipids.itp"))
    print("Moleculetypes:", [m.name for m in mols])
    for m in mols:
        print(m.bonds)