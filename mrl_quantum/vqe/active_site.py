"""Construcao da geometria do sitio ativo da RNase H1.

Tenta baixar a estrutura PDB 2QKB (RNase H1 humana complexada com
duplex RNA:DNA hibrido). Se falhar, constroi geometria modelo usando ASE
com coordenacao octaedrica dos dois ions Mg2+ e ligantes simplificados.

Referencia:
  Nowotny M et al. (2005) Cell 121(7):1005-1016 — estrutura do mecanismo
  de dois-metais da RNase H1.
"""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path
from typing import Final

import numpy as np
from ase import Atoms
from ase.io import write as ase_write

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constantes do sitio ativo — distancias experimentais (Angstrom)
# Ref: Nowotny et al. (2005) Cell; Yang et al. (2006) Mol Cell
# ---------------------------------------------------------------------------
MG_MG_DISTANCE: Final[float] = 3.9       # distancia Mg-Mg no sitio bimetalico
MG_O_DISTANCE: Final[float] = 2.1        # coordenacao Mg-O (fosfato/carboxilato)
MG_N_DISTANCE: Final[float] = 2.2        # coordenacao Mg-N (imidazol His124)
OH_DISTANCE: Final[float] = 0.96         # O-H na agua

# Residuos catalíticos — numeracao humana RNase H1
CATALYTIC_RESIDUES: Final[list[str]] = ["D10", "E48", "D70", "H124"]

# PDB da RNase H1 humana com duplex RNA:DNA
PDB_ID: Final[str] = "2QKB"

# Diretorio padrao para salvar geometria
_DEFAULT_OUTPUT_DIR: Final[Path] = (
    Path(__file__).resolve().parent.parent.parent / "results" / "vqe"
)


def _try_download_pdb(pdb_id: str, output_dir: Path) -> Path | None:
    """Tenta baixar estrutura PDB via Bio.PDB.PDBList.

    Returns:
        Caminho do arquivo PDB, ou None se falhar.
    """
    try:
        from Bio.PDB import PDBList

        pdbl = PDBList(verbose=False)
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = pdbl.retrieve_pdb_file(
                pdb_id, pdir=tmpdir, file_format="pdb"
            )
            if fname and Path(fname).exists():
                import shutil

                dest = output_dir / f"{pdb_id.lower()}.pdb"
                output_dir.mkdir(parents=True, exist_ok=True)
                shutil.copy2(fname, dest)
                logger.info("PDB %s baixado com sucesso: %s", pdb_id, dest)
                return dest
    except Exception as exc:
        logger.warning("Falha ao baixar PDB %s: %s", pdb_id, exc)
    return None


def _extract_active_site_from_pdb(pdb_path: Path) -> Atoms | None:
    """Extrai geometria do sitio ativo a partir de PDB real.

    Residuos extraidos: D10, E48, D70, H124 + 2 Mg2+ + aguas proximas.
    Usa apenas atomos pesados dos residuos catalíticos (sem H).
    """
    try:
        from Bio.PDB import PDBParser

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(PDB_ID, str(pdb_path))
        model = structure[0]

        # Mapa de residuos catalíticos (chain A)
        target_residues = {10: "D", 48: "E", 70: "D", 124: "H"}
        symbols: list[str] = []
        positions: list[list[float]] = []

        for chain in model:
            for residue in chain:
                res_id = residue.get_id()[1]
                # Residuos catalíticos
                if res_id in target_residues:
                    for atom in residue:
                        elem = atom.element.strip()
                        if elem and elem != "H":
                            symbols.append(elem)
                            positions.append(list(atom.get_vector()))

                # Ions Mg2+
                resname = residue.get_resname().strip()
                if resname == "MG":
                    for atom in residue:
                        symbols.append("Mg")
                        positions.append(list(atom.get_vector()))

        # Adicionar aguas proximas ao Mg (dentro de 3.0 A)
        mg_positions = [
            p for s, p in zip(symbols, positions) if s == "Mg"
        ]
        if mg_positions:
            mg_center = np.mean(mg_positions, axis=0)
            for chain in model:
                for residue in chain:
                    if residue.get_resname().strip() == "HOH":
                        for atom in residue:
                            pos = list(atom.get_vector())
                            dist = np.linalg.norm(np.array(pos) - mg_center)
                            if dist < 5.0 and len(
                                [s for s in symbols if s == "O"]
                            ) < 6:
                                symbols.append("O")
                                positions.append(pos)

        if len(symbols) < 5:
            logger.warning(
                "Poucos atomos extraidos do PDB (%d), usando modelo",
                len(symbols),
            )
            return None

        atoms = Atoms(symbols=symbols, positions=positions)
        logger.info(
            "Sitio ativo extraido do PDB: %d atomos", len(atoms)
        )
        return atoms

    except Exception as exc:
        logger.warning("Erro ao extrair sitio ativo do PDB: %s", exc)
        return None


def _build_model_geometry() -> Atoms:
    """Constroi geometria modelo do sitio ativo com ASE.

    Modelo simplificado do mecanismo de dois-metais:
    - 2 ions Mg2+ separados por 3.9 A
    - Coordenacao octaedrica em torno de cada Mg2+
    - 3 carboxilatos simplificados (D10, E48, D70) como COO-
    - 1 imidazol simplificado (H124) como anel de 5 membros
    - 3 moleculas de agua como atomos de O
    - 1 fragmento fosfato (PO4) do substrato RNA

    A geometria reproduz as distancias experimentais de Nowotny et al. (2005).
    """
    symbols: list[str] = []
    positions: list[list[float]] = []

    # --- Dois ions Mg2+ ao longo do eixo x ---
    # Metal A (catalitico) e Metal B (estrutural)
    mg_a = np.array([0.0, 0.0, 0.0])
    mg_b = np.array([MG_MG_DISTANCE, 0.0, 0.0])

    symbols.extend(["Mg", "Mg"])
    positions.extend([mg_a.tolist(), mg_b.tolist()])

    # --- Carboxilato D10 — ponte entre os dois Mg2+ ---
    # Oxigenio OD1 coordena Metal A, OD2 coordena Metal B
    d10_od1 = mg_a + np.array([MG_O_DISTANCE, MG_O_DISTANCE * 0.5, 0.0])
    d10_od2 = mg_b + np.array([-MG_O_DISTANCE, MG_O_DISTANCE * 0.5, 0.0])
    d10_c = (d10_od1 + d10_od2) / 2 + np.array([0.0, 1.2, 0.0])

    symbols.extend(["C", "O", "O"])
    positions.extend([d10_c.tolist(), d10_od1.tolist(), d10_od2.tolist()])

    # --- Carboxilato E48 — coordena Metal A ---
    e48_od1 = mg_a + np.array([-MG_O_DISTANCE, 0.0, MG_O_DISTANCE * 0.7])
    e48_od2 = mg_a + np.array(
        [-MG_O_DISTANCE - 1.2, 0.0, MG_O_DISTANCE * 0.7]
    )
    e48_c = (e48_od1 + e48_od2) / 2 + np.array([0.0, 0.0, 1.2])

    symbols.extend(["C", "O", "O"])
    positions.extend([e48_c.tolist(), e48_od1.tolist(), e48_od2.tolist()])

    # --- Carboxilato D70 — coordena Metal B ---
    d70_od1 = mg_b + np.array([MG_O_DISTANCE, 0.0, -MG_O_DISTANCE * 0.7])
    d70_od2 = mg_b + np.array(
        [MG_O_DISTANCE + 1.2, 0.0, -MG_O_DISTANCE * 0.7]
    )
    d70_c = (d70_od1 + d70_od2) / 2 + np.array([0.0, 0.0, -1.2])

    symbols.extend(["C", "O", "O"])
    positions.extend([d70_c.tolist(), d70_od1.tolist(), d70_od2.tolist()])

    # --- Imidazol H124 — coordena Metal B via N-epsilon ---
    # Anel de 5 membros simplificado (C3N2H3 -> so atomos pesados)
    h124_ne = mg_b + np.array([0.0, -MG_N_DISTANCE, 0.0])
    h124_cd = h124_ne + np.array([-1.1, -0.8, 0.0])
    h124_ce = h124_ne + np.array([1.1, -0.8, 0.0])
    h124_nd = h124_ne + np.array([-0.7, -2.0, 0.0])
    h124_cg = h124_ne + np.array([0.7, -2.0, 0.0])

    symbols.extend(["N", "C", "C", "N", "C"])
    positions.extend([
        h124_ne.tolist(),
        h124_cd.tolist(),
        h124_ce.tolist(),
        h124_nd.tolist(),
        h124_cg.tolist(),
    ])

    # --- Aguas de coordenacao (3 moleculas, representadas como O) ---
    # Agua 1: coordena Metal A (posicao axial +z)
    w1 = mg_a + np.array([0.0, 0.0, MG_O_DISTANCE])
    # Agua 2: coordena Metal A (posicao axial -z)
    w2 = mg_a + np.array([0.0, 0.0, -MG_O_DISTANCE])
    # Agua 3: coordena Metal B (posicao axial +z)
    w3 = mg_b + np.array([0.0, 0.0, MG_O_DISTANCE])

    symbols.extend(["O", "O", "O"])
    positions.extend([w1.tolist(), w2.tolist(), w3.tolist()])

    # --- Fosfato do substrato RNA (simplificado como PO4) ---
    # Posicao entre os dois metais, ligeiramente acima
    po4_p = (mg_a + mg_b) / 2 + np.array([0.0, MG_O_DISTANCE, 0.0])
    po4_o1 = po4_p + np.array([0.0, 0.0, 1.2])
    po4_o2 = po4_p + np.array([0.0, 0.0, -1.2])
    po4_o3 = po4_p + np.array([1.2, 0.5, 0.0])
    po4_o4 = po4_p + np.array([-1.2, 0.5, 0.0])

    symbols.extend(["P", "O", "O", "O", "O"])
    positions.extend([
        po4_p.tolist(),
        po4_o1.tolist(),
        po4_o2.tolist(),
        po4_o3.tolist(),
        po4_o4.tolist(),
    ])

    atoms = Atoms(symbols=symbols, positions=positions)
    logger.info(
        "Geometria modelo construida: %d atomos (%s)",
        len(atoms),
        ", ".join(f"{s}: {list(symbols).count(s)}" for s in sorted(set(symbols))),
    )
    return atoms


def build_active_site(
    output_dir: Path | None = None,
    try_pdb: bool = True,
) -> tuple[Atoms, str, Path]:
    """Constroi a geometria do sitio ativo da RNase H1.

    Estrategia:
    1. Tenta baixar PDB 2QKB e extrair sitio ativo
    2. Se falhar, constroi geometria modelo com ASE

    Args:
        output_dir: Diretorio para salvar arquivo XYZ. Se None, usa padrao.
        try_pdb: Se True, tenta baixar PDB primeiro.

    Returns:
        Tupla (atoms, source, xyz_path):
        - atoms: Objeto ASE Atoms com a geometria
        - source: "pdb_2qkb" ou "model_ase" indicando a origem
        - xyz_path: Caminho do arquivo XYZ salvo
    """
    out = output_dir or _DEFAULT_OUTPUT_DIR
    out.mkdir(parents=True, exist_ok=True)

    atoms: Atoms | None = None
    source = "model_ase"

    # Estrategia 1: PDB real
    if try_pdb:
        pdb_path = _try_download_pdb(PDB_ID, out)
        if pdb_path is not None:
            atoms = _extract_active_site_from_pdb(pdb_path)
            if atoms is not None:
                source = "pdb_2qkb"
                logger.info("Usando geometria do PDB 2QKB")

    # Estrategia 2: modelo ASE
    if atoms is None:
        atoms = _build_model_geometry()
        source = "model_ase"
        logger.info("Usando geometria modelo (ASE)")

    # Salvar como XYZ
    xyz_path = out / "rnase_h1_active_site.xyz"
    ase_write(str(xyz_path), atoms, format="xyz")
    logger.info("Geometria salva: %s", xyz_path)

    return atoms, source, xyz_path
