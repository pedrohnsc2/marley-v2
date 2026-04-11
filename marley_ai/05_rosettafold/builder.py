"""Construcao de modelos 3D de duplex RNA:DNA hibrido (ASO:SL RNA).

Gera coordenadas atomicas para um duplex idealizado na forma-A,
aplicando modificacoes LNA e backbone fosforotioato (PS).

O modelo gerado pode ser aberto em PyMOL/ChimeraX para visualizacao.

Ref: Saenger W (1984) Principles of Nucleic Acid Structure.
     Egli M et al. (2005) Chem Biol 12:669-675 (LNA constraints).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Final

import importlib as _il

import numpy as np

# Importacao via importlib — Python nao permite "from marley_ai.05_rosettafold..."
# porque "05" e interpretado como literal numerico invalido no parser.
_cfg = _il.import_module("marley_ai.05_rosettafold.config")

A_FORM_INCLINATION = _cfg.A_FORM_INCLINATION
A_FORM_RISE = _cfg.A_FORM_RISE
A_FORM_TWIST = _cfg.A_FORM_TWIST
A_FORM_X_DISPLACEMENT = _cfg.A_FORM_X_DISPLACEMENT
ASO_SEQUENCE = _cfg.ASO_SEQUENCE
COMPLEMENT = _cfg.COMPLEMENT
DNA_SUGAR_PUCKER_PHASE = _cfg.DNA_SUGAR_PUCKER_PHASE
LNA_INFLUENCE_RADIUS = _cfg.LNA_INFLUENCE_RADIUS
LNA_POSITIONS = _cfg.LNA_POSITIONS
LNA_SUGAR_PUCKER_PHASE = _cfg.LNA_SUGAR_PUCKER_PHASE
PS_BOND_ANGLE_SHIFT = _cfg.PS_BOND_ANGLE_SHIFT
PS_BOND_LENGTH = _cfg.PS_BOND_LENGTH
PO_BOND_LENGTH = _cfg.PO_BOND_LENGTH
RNA_COMPLEMENT = _cfg.RNA_COMPLEMENT
SL_RNA_SEQUENCE = _cfg.SL_RNA_SEQUENCE
TARGET_RNA_SEQUENCE = _cfg.TARGET_RNA_SEQUENCE
ASO_TARGET_START = _cfg.ASO_TARGET_START
ASO_TARGET_END = _cfg.ASO_TARGET_END
StructuralAnalysisConfig = _cfg.StructuralAnalysisConfig


# ---------------------------------------------------------------------------
# Coordenadas atomicas idealizadas de bases nitrogenadas
# Sistema de referencia: base no plano XY, ponto de attache (N1 ou N9) na origem
# Ref: Parkinson G et al. (1996) Acta Cryst D52:57-64
# ---------------------------------------------------------------------------

# Purinas (A, G) — ponto de attach no N9
_PURINE_ATOMS: Final[dict[str, dict[str, tuple[float, float, float]]]] = {
    "A": {
        "N9":  ( 0.000,  0.000, 0.0),
        "C8":  ( 1.200,  0.700, 0.0),
        "N7":  ( 2.100,  0.000, 0.0),
        "C5":  ( 1.500, -1.200, 0.0),
        "C6":  ( 2.000, -2.400, 0.0),
        "N6":  ( 3.300, -2.600, 0.0),
        "N1":  ( 1.200, -3.500, 0.0),
        "C2":  ( 0.000, -3.300, 0.0),
        "N3":  (-0.500, -2.100, 0.0),
        "C4":  ( 0.300, -1.000, 0.0),
    },
    "G": {
        "N9":  ( 0.000,  0.000, 0.0),
        "C8":  (-1.037,  0.925, 0.0),
        "N7":  (-0.649,  1.998, 0.0),
        "C5":  ( 0.678,  1.798, 0.0),
        "C6":  ( 1.664,  2.644, 0.0),
        "O6":  ( 1.484,  3.847, 0.0),
        "N1":  ( 2.958,  2.223, 0.0),
        "C2":  ( 3.138,  1.020, 0.0),
        "N2":  ( 4.306,  0.663, 0.0),
        "N3":  ( 2.152,  0.174, 0.0),
        "C4":  ( 0.858,  0.595, 0.0),
    },
}

# Pirimidinas (C, T/U) — ponto de attach no N1
_PYRIMIDINE_ATOMS: Final[dict[str, dict[str, tuple[float, float, float]]]] = {
    "C": {
        "N1":  ( 0.000,  0.000, 0.0),
        "C2":  ( 0.677,  1.110, 0.0),
        "O2":  ( 0.214,  2.257, 0.0),
        "N3":  ( 1.993,  0.830, 0.0),
        "C4":  ( 2.514, -0.398, 0.0),
        "N4":  ( 3.610, -0.715, 0.0),
        "C5":  ( 1.698, -1.486, 0.0),
        "C6":  ( 0.463, -1.147, 0.0),
    },
    "T": {
        "N1":  ( 0.000,  0.000, 0.0),
        "C2":  (-1.265,  0.301, 0.0),
        "O2":  (-2.213, -0.493, 0.0),
        "N3":  (-1.405,  1.639, 0.0),
        "C4":  (-0.398,  2.515, 0.0),
        "O4":  (-0.436,  3.654, 0.0),
        "C5":  ( 0.888,  2.074, 0.0),
        "C7":  ( 1.976,  2.891, 0.0),  # grupo metila da timina
        "C6":  ( 0.948,  0.795, 0.0),
    },
    "U": {
        "N1":  ( 0.000,  0.000, 0.0),
        "C2":  (-1.265,  0.301, 0.0),
        "O2":  (-2.213, -0.493, 0.0),
        "N3":  (-1.405,  1.639, 0.0),
        "C4":  (-0.398,  2.515, 0.0),
        "O4":  (-0.436,  3.654, 0.0),
        "C5":  ( 0.888,  2.074, 0.0),
        "C6":  ( 0.948,  0.795, 0.0),
    },
}

# Unificar em um unico dicionario
BASE_ATOMS: Final[dict[str, dict[str, tuple[float, float, float]]]] = {
    **_PURINE_ATOMS,
    **_PYRIMIDINE_ATOMS,
}

# Classificacao: purina ou pirimidina
PURINES: Final[set[str]] = {"A", "G"}
PYRIMIDINES: Final[set[str]] = {"C", "T", "U"}


# ---------------------------------------------------------------------------
# Estrutura de dados para coordenadas
# ---------------------------------------------------------------------------

@dataclass
class Atom:
    """Representa um atomo em coordenadas cartesianas."""
    name: str           # nome do atomo (ex: "N9", "P", "C1'")
    residue_name: str   # nome do residuo (ex: "A", "DA", "RC")
    chain_id: str       # cadeia (A = fita sense, B = fita antisense)
    residue_num: int    # numero do residuo (1-indexed)
    x: float
    y: float
    z: float
    element: str        # simbolo do elemento (C, N, O, P, S)
    b_factor: float = 0.0  # fator B — usado para indicar LNA/PS


# ---------------------------------------------------------------------------
# Funcoes auxiliares de geometria
# ---------------------------------------------------------------------------

def _rotation_matrix_z(theta_deg: float) -> np.ndarray:
    """Matriz de rotacao ao redor do eixo Z.

    Args:
        theta_deg: Angulo de rotacao em graus.

    Returns:
        Matriz 3x3 de rotacao.
    """
    theta = math.radians(theta_deg)
    cos_t = math.cos(theta)
    sin_t = math.sin(theta)
    return np.array([
        [cos_t, -sin_t, 0.0],
        [sin_t,  cos_t, 0.0],
        [0.0,    0.0,   1.0],
    ])


def _rotation_matrix_y(theta_deg: float) -> np.ndarray:
    """Matriz de rotacao ao redor do eixo Y.

    Args:
        theta_deg: Angulo de rotacao em graus.

    Returns:
        Matriz 3x3 de rotacao.
    """
    theta = math.radians(theta_deg)
    cos_t = math.cos(theta)
    sin_t = math.sin(theta)
    return np.array([
        [ cos_t, 0.0, sin_t],
        [ 0.0,   1.0, 0.0  ],
        [-sin_t, 0.0, cos_t],
    ])


def _rotation_matrix_x(theta_deg: float) -> np.ndarray:
    """Matriz de rotacao ao redor do eixo X.

    Args:
        theta_deg: Angulo de rotacao em graus.

    Returns:
        Matriz 3x3 de rotacao.
    """
    theta = math.radians(theta_deg)
    cos_t = math.cos(theta)
    sin_t = math.sin(theta)
    return np.array([
        [1.0, 0.0,    0.0   ],
        [0.0, cos_t, -sin_t ],
        [0.0, sin_t,  cos_t ],
    ])


def _get_residue_name(base: str, chain_type: str, is_lna: bool = False) -> str:
    """Gera nome de residuo PDB para uma base.

    Convencao:
        - DNA: DA, DC, DG, DT
        - RNA: A, C, G, U
        - LNA: marcado com B-factor alto (nao tem codigo PDB especifico)

    Args:
        base: Letra da base (A, C, G, T, U).
        chain_type: "DNA" ou "RNA".
        is_lna: Se True, marca como modificacao LNA.

    Returns:
        Nome do residuo em formato PDB.
    """
    if chain_type == "DNA":
        # Nome PDB para desoxirribonucleotideo
        return f"D{base}" if base != "U" else "DT"
    # RNA — nome de 1 ou 2 caracteres
    if base == "T":
        return "U"
    return base


# ---------------------------------------------------------------------------
# Construtor principal de duplex
# ---------------------------------------------------------------------------

def build_duplex(
    rna_sequence: str,
    dna_sequence: str,
    *,
    rise: float = A_FORM_RISE,
    twist: float = A_FORM_TWIST,
    inclination: float = A_FORM_INCLINATION,
    x_displacement: float = A_FORM_X_DISPLACEMENT,
    lna_positions: list[int] | None = None,
    use_ps_backbone: bool = True,
) -> list[Atom]:
    """Constroi duplex RNA:DNA hibrido idealizado na forma-A.

    A fita RNA (SL RNA) e a cadeia A, o ASO (DNA) e a cadeia B.
    As bases sao posicionadas usando transformacoes helicoidais:
    cada par de bases e rotacionado (twist) e transladado (rise) ao
    longo do eixo Z, com inclinacao (inclination) e deslocamento-x.

    Geometria do par de bases:
    - O C1' da fita RNA esta no raio positivo (backbone externo)
    - O C1' da fita DNA esta no raio oposto (backbone externo)
    - As bases se estendem para dentro, em direcao ao centro da helice
    - Distancia N-N entre bases pareadas: ~2.8 A (H-bond distance)
    - Distancia C1'-C1': ~10.4 A

    Args:
        rna_sequence: Sequencia da fita RNA (5'->3', letras ACGU).
        dna_sequence: Sequencia da fita DNA/ASO (5'->3', letras ACGT).
        rise: Ascensao por par de bases (Angstroms).
        twist: Rotacao por par de bases (graus).
        inclination: Inclinacao das bases (graus).
        x_displacement: Deslocamento-x das bases (Angstroms).
        lna_positions: Posicoes LNA no ASO (0-indexed). None = sem LNA.
        use_ps_backbone: Se True, aplica distorcao PS no backbone do ASO.

    Returns:
        Lista de Atom com todas as coordenadas do duplex.
    """
    if lna_positions is None:
        lna_positions = []

    atoms: list[Atom] = []
    n_bp = min(len(rna_sequence), len(dna_sequence))

    # Distancia C1'-C1' entre fitas opostas (Angstroms)
    c1_c1_distance = 10.4

    # Raio do backbone (distancia do eixo helicoidal ao C1')
    backbone_radius = c1_c1_distance / 2.0  # 5.2 A

    # Distancia do C1' ao ponto de attach da base (N9 purina, N1 pirimidina)
    # Comprimento da ligacao glicosidica C1'-N9/N1: ~1.47 A
    glycosidic_length = 1.47

    # Escala das coordenadas da base — 1.0 para manter geometria realista
    # Com scale=1.0, as bases de Watson-Crick ficam a ~2.8 A (distancia H-bond)
    base_scale = 1.0

    for i in range(n_bp):
        rna_base = rna_sequence[i]
        dna_base = dna_sequence[i]

        # --- Transformacao helicoidal ---
        total_twist = twist * i
        z_offset = rise * i

        rot_z = _rotation_matrix_z(total_twist)
        rot_incl = _rotation_matrix_x(inclination)
        rot_total = rot_z @ rot_incl

        # Centro do par de bases (com x_displacement da forma-A)
        center = rot_z @ np.array([x_displacement, 0.0, 0.0])
        center[2] += z_offset

        # --- Fita RNA (cadeia A) ---
        # C1' no lado positivo do raio (backbone externo)
        c1_rna_local = np.array([backbone_radius, 0.0, 0.0])
        c1_rna = center + rot_z @ c1_rna_local

        is_purine_rna = rna_base in PURINES
        rna_base_key = rna_base if rna_base in BASE_ATOMS else "A"
        rna_res_name = _get_residue_name(rna_base, "RNA")

        # Fosfato P: deslocado ao longo do backbone (mais externo que C1')
        if i > 0:
            p_rna_local = np.array([backbone_radius + 1.6, 2.4, 0.0])
            p_pos_rna = center + rot_z @ p_rna_local
            atoms.append(Atom(
                name="P", residue_name=rna_res_name, chain_id="A",
                residue_num=i + 1,
                x=float(p_pos_rna[0]), y=float(p_pos_rna[1]),
                z=float(p_pos_rna[2]),
                element="P",
            ))

        # C4' (entre C1' e P, levemente externo)
        c4_rna_local = np.array([backbone_radius + 0.8, 1.2, 0.0])
        c4_pos_rna = center + rot_z @ c4_rna_local
        atoms.append(Atom(
            name="C4'", residue_name=rna_res_name, chain_id="A",
            residue_num=i + 1,
            x=float(c4_pos_rna[0]), y=float(c4_pos_rna[1]),
            z=float(c4_pos_rna[2]),
            element="C",
        ))

        # Atomos da base RNA — orientados PARA DENTRO (em direcao ao centro)
        # O attach point (N9/N1) esta em C1' - glycosidic_length na direcao do centro
        base_atoms_dict = BASE_ATOMS.get(rna_base_key, BASE_ATOMS["A"])
        for atom_name, local_coords in base_atoms_dict.items():
            # Coordenadas locais da base: escalar e orientar para dentro
            local = np.array(local_coords) * base_scale
            # Bases RNA apontam para -X (para dentro da helice)
            local[0] = -abs(local[0])
            # Posicionar a partir do C1', deslocado para dentro
            world = c1_rna - rot_z @ np.array([glycosidic_length, 0.0, 0.0])
            world = world + rot_total @ local
            elem = atom_name[0]
            atoms.append(Atom(
                name=atom_name, residue_name=rna_res_name, chain_id="A",
                residue_num=i + 1,
                x=float(world[0]), y=float(world[1]), z=float(world[2]),
                element=elem,
            ))

        # --- Fita DNA/ASO (cadeia B) ---
        # C1' no lado negativo do raio (backbone oposto)
        c1_dna_local = np.array([-backbone_radius, 0.0, 0.0])
        c1_dna = center + rot_z @ c1_dna_local

        is_lna = i in lna_positions
        is_purine_dna = dna_base in PURINES
        dna_base_key = dna_base if dna_base in BASE_ATOMS else "T"
        dna_res_name = _get_residue_name(dna_base, "DNA", is_lna=is_lna)

        # B-factor: 20.0 para LNA, 10.0 para PS, 0.0 para normal
        b_factor = 20.0 if is_lna else (10.0 if use_ps_backbone else 0.0)

        # Fosfato P da fita DNA
        ps_shift = PS_BOND_ANGLE_SHIFT if use_ps_backbone else 0.0
        if i > 0:
            p_dna_local = np.array([
                -(backbone_radius + 1.6),
                -2.4 + math.radians(ps_shift),
                0.0,
            ])
            p_pos_dna = center + rot_z @ p_dna_local
            atoms.append(Atom(
                name="P", residue_name=dna_res_name, chain_id="B",
                residue_num=i + 1,
                x=float(p_pos_dna[0]), y=float(p_pos_dna[1]),
                z=float(p_pos_dna[2]),
                element="P", b_factor=b_factor,
            ))

        # C4' do acucar DNA
        sugar_phase = LNA_SUGAR_PUCKER_PHASE if is_lna else DNA_SUGAR_PUCKER_PHASE
        sugar_displacement = 0.3 * math.cos(math.radians(sugar_phase))
        c4_dna_local = np.array([
            -(backbone_radius + 0.8),
            -1.2 + sugar_displacement,
            0.0,
        ])
        c4_pos_dna = center + rot_z @ c4_dna_local
        atoms.append(Atom(
            name="C4'", residue_name=dna_res_name, chain_id="B",
            residue_num=i + 1,
            x=float(c4_pos_dna[0]), y=float(c4_pos_dna[1]),
            z=float(c4_pos_dna[2]),
            element="C", b_factor=b_factor,
        ))

        # Atomos da base DNA — orientados PARA DENTRO (em direcao a base RNA)
        base_atoms_dict_dna = BASE_ATOMS.get(dna_base_key, BASE_ATOMS["T"])
        for atom_name, local_coords in base_atoms_dict_dna.items():
            local = np.array(local_coords) * base_scale
            # Bases DNA apontam para +X (para dentro, em direcao a fita RNA)
            local[0] = abs(local[0])
            # Posicionar a partir do C1' DNA, deslocado para dentro (+X)
            world = c1_dna + rot_z @ np.array([glycosidic_length, 0.0, 0.0])
            world = world + rot_total @ local
            elem = atom_name[0]
            atoms.append(Atom(
                name=atom_name, residue_name=dna_res_name, chain_id="B",
                residue_num=i + 1,
                x=float(world[0]), y=float(world[1]), z=float(world[2]),
                element=elem, b_factor=b_factor,
            ))

    return atoms


def build_single_strand_rna(
    sequence: str,
    *,
    rise: float = A_FORM_RISE,
    twist: float = A_FORM_TWIST,
    chain_id: str = "A",
) -> list[Atom]:
    """Constroi modelo helicoidal de uma fita simples de RNA.

    Usado para modelar o SL RNA isolado (sem ASO ligado).
    A geometria segue forma-A mas sem a fita complementar.

    Args:
        sequence: Sequencia RNA (5'->3', letras ACGU).
        rise: Ascensao por residuo (Angstroms).
        twist: Rotacao por residuo (graus).
        chain_id: Identificador da cadeia PDB.

    Returns:
        Lista de Atom com coordenadas do RNA de fita simples.
    """
    atoms: list[Atom] = []
    radius = 10.0  # raio da helice para fita simples

    for i, base in enumerate(sequence):
        total_twist = twist * i
        z_offset = rise * i
        rot_z = _rotation_matrix_z(total_twist)

        # Posicao ao longo da helice
        origin = rot_z @ np.array([radius, 0.0, 0.0])
        origin[2] += z_offset

        res_name = _get_residue_name(base, "RNA")
        is_purine = base in PURINES

        # Backbone P
        if i > 0:
            p_offset = rot_z @ np.array([radius + 1.6, 2.4, 0.0])
            p_offset[2] += z_offset
            atoms.append(Atom(
                name="P", residue_name=res_name, chain_id=chain_id,
                residue_num=i + 1,
                x=float(p_offset[0]), y=float(p_offset[1]), z=float(p_offset[2]),
                element="P",
            ))

        # C4'
        c4_offset = rot_z @ np.array([radius + 0.8, 1.2, 0.0])
        c4_offset[2] += z_offset
        atoms.append(Atom(
            name="C4'", residue_name=res_name, chain_id=chain_id,
            residue_num=i + 1,
            x=float(c4_offset[0]), y=float(c4_offset[1]), z=float(c4_offset[2]),
            element="C",
        ))

        # Atomos da base
        base_key = base if base in BASE_ATOMS else "A"
        base_atoms_dict = BASE_ATOMS[base_key]
        rot_incl = _rotation_matrix_x(A_FORM_INCLINATION)
        rot_total = rot_z @ rot_incl

        for atom_name, local_coords in base_atoms_dict.items():
            local = np.array(local_coords)
            world = origin + rot_total @ local
            elem = atom_name[0]
            atoms.append(Atom(
                name=atom_name, residue_name=res_name, chain_id=chain_id,
                residue_num=i + 1,
                x=float(world[0]), y=float(world[1]), z=float(world[2]),
                element=elem,
            ))

    return atoms


# ---------------------------------------------------------------------------
# Exportacao PDB
# ---------------------------------------------------------------------------

def atoms_to_pdb(
    atoms: list[Atom],
    *,
    title: str = "ASO:SL RNA Duplex Model",
    remarks: list[str] | None = None,
) -> str:
    """Converte lista de atomos em texto PDB valido.

    Segue o formato PDB padrao (colunas fixas):
    ATOM  serial name resName chainID resSeq x y z occupancy bFactor element

    Ref: https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html

    Args:
        atoms: Lista de Atom com coordenadas.
        title: Titulo para o cabecalho TITLE.
        remarks: Linhas REMARK adicionais.

    Returns:
        String com conteudo PDB completo.
    """
    lines: list[str] = []

    # Cabecalho
    lines.append(f"TITLE     {title}")
    lines.append(
        f"REMARK   1 Generated by Marley - Structural Analysis Module"
    )
    lines.append(
        f"REMARK   2 Date: {datetime.now(tz=timezone.utc).isoformat()}"
    )
    if remarks:
        for i, remark in enumerate(remarks, start=3):
            lines.append(f"REMARK  {i:2d} {remark}")

    # Registros ATOM
    # Formato PDB: colunas fixas
    # 1-6: "ATOM  "
    # 7-11: serial (int, right-justified)
    # 13-16: name (atom name)
    # 17: altLoc
    # 18-20: resName
    # 22: chainID
    # 23-26: resSeq (int, right-justified)
    # 31-38: x (float 8.3)
    # 39-46: y (float 8.3)
    # 47-54: z (float 8.3)
    # 55-60: occupancy (float 6.2)
    # 61-66: tempFactor (float 6.2)
    # 77-78: element

    # Ordenar atomos por cadeia e numero de residuo para PDB valido
    # (cadeias devem ser contiguas no formato PDB)
    sorted_atoms = sorted(atoms, key=lambda a: (a.chain_id, a.residue_num))

    serial = 1
    prev_chain = ""
    for atom in sorted_atoms:
        # Inserir TER entre cadeias
        if prev_chain and atom.chain_id != prev_chain:
            lines.append("TER")
        prev_chain = atom.chain_id

        # Formatacao do nome do atomo: 4 caracteres
        # Atomos com 1 caractere de elemento comecam na coluna 14
        # Atomos com 2 caracteres comecam na coluna 13
        if len(atom.name) <= 3:
            atom_name_fmt = f" {atom.name:<3s}"
        else:
            atom_name_fmt = f"{atom.name:<4s}"

        # Formatacao do nome do residuo: 3 caracteres, right-justified
        res_name_fmt = f"{atom.residue_name:>3s}"

        line = (
            f"ATOM  {serial:5d} {atom_name_fmt}"
            f" {res_name_fmt} {atom.chain_id}{atom.residue_num:4d}    "
            f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
            f"{1.00:6.2f}{atom.b_factor:6.2f}"
            f"          {atom.element:>2s}  "
        )
        lines.append(line)
        serial += 1

    lines.append("TER")
    lines.append("END")

    return "\n".join(lines)


def save_pdb(
    atoms: list[Atom],
    output_path: Path,
    *,
    title: str = "ASO:SL RNA Duplex Model",
    remarks: list[str] | None = None,
) -> Path:
    """Salva lista de atomos como arquivo PDB.

    Args:
        atoms: Lista de Atom com coordenadas.
        output_path: Caminho do arquivo de saida.
        title: Titulo para o cabecalho TITLE.
        remarks: Linhas REMARK adicionais.

    Returns:
        Caminho do arquivo salvo.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_text = atoms_to_pdb(atoms, title=title, remarks=remarks)
    output_path.write_text(pdb_text, encoding="utf-8")
    return output_path


# ---------------------------------------------------------------------------
# Funcao principal de construcao
# ---------------------------------------------------------------------------

def build_aso_sl_complex(
    config: StructuralAnalysisConfig | None = None,
) -> dict[str, list[Atom]]:
    """Constroi modelos 3D completos para analise estrutural.

    Gera tres modelos:
    1. SL RNA isolado (fita simples helicoidal, 39 nt)
    2. Duplex ASO:SL RNA (regiao de hibridizacao, 24 bp)
    3. SL RNA com ASO ligado (39 nt RNA + 24 nt ASO na regiao-alvo)

    Args:
        config: Configuracao do modulo. Se None, usa defaults.

    Returns:
        Dicionario com tres conjuntos de atomos:
        - "sl_rna_free": SL RNA isolado
        - "duplex": Duplex ASO:RNA na regiao de hibridizacao
        - "complex": SL RNA completo com ASO ligado
    """
    if config is None:
        config = StructuralAnalysisConfig()

    # Extrair a regiao do SL RNA que hibridiza com o ASO
    # O ASO e antisense (5'->3') em relacao a regiao-alvo do SL RNA
    target_rna = SL_RNA_SEQUENCE[ASO_TARGET_START:ASO_TARGET_END]

    # O ASO hibrida com a regiao complementar do SL RNA
    # RNA sequence: 5'-AACGCUAUAUAAGUAUCAGUUUCUG-3' (posicoes 5-29 do SL RNA, 25 nt)
    # ASO sequence: 3'-TGCGATATATCATAGTTAAAGAC...-5' (leitura 3'->5')
    # O ASO como lido 5'->3' e a sequencia que temos: ACAGAAACTGATACTTATATAGCGT

    # Para o duplex, alinhamos RNA 5'->3' com ASO 3'->5'
    # Ou seja, RNA[i] pareia com ASO[n-1-i]
    # Vamos construir o duplex na orientacao padrao:
    # RNA (cadeia A) de 5' para 3', ASO (cadeia B) de 3' para 5'
    aso_reversed = ASO_SEQUENCE[::-1]  # reverter para alinhar com RNA

    # Construir SL RNA isolado
    sl_rna_free = build_single_strand_rna(SL_RNA_SEQUENCE)

    # Construir duplex hibridizado
    # Mapear posicoes LNA do ASO original para o ASO revertido
    aso_len = len(ASO_SEQUENCE)
    lna_reversed = [aso_len - 1 - pos for pos in LNA_POSITIONS]

    duplex = build_duplex(
        rna_sequence=target_rna,
        dna_sequence=aso_reversed,
        rise=config.rise_per_bp,
        twist=config.twist_per_bp,
        inclination=config.inclination,
        x_displacement=config.x_displacement,
        lna_positions=lna_reversed,
        use_ps_backbone=config.use_ps_backbone,
    )

    # Construir complexo completo: SL RNA + ASO na regiao-alvo
    # O SL RNA e de fita simples nos flancos, duplex na regiao-alvo
    complex_atoms: list[Atom] = []

    # Flanco 5' do SL RNA (posicoes 0-4, fita simples)
    flank_5 = build_single_strand_rna(
        SL_RNA_SEQUENCE[:ASO_TARGET_START],
        chain_id="A",
    )
    complex_atoms.extend(flank_5)

    # Regiao de hibridizacao — re-numerar para continuidade
    max_resnum_a = ASO_TARGET_START  # ultimo residuo do flanco 5'
    for atom in duplex:
        new_atom = Atom(
            name=atom.name,
            residue_name=atom.residue_name,
            chain_id=atom.chain_id,
            residue_num=(atom.residue_num + max_resnum_a
                         if atom.chain_id == "A"
                         else atom.residue_num),
            x=atom.x,
            y=atom.y,
            z=atom.z + ASO_TARGET_START * config.rise_per_bp,
            element=atom.element,
            b_factor=atom.b_factor,
        )
        complex_atoms.append(new_atom)

    # Flanco 3' do SL RNA (posicoes 30-38, fita simples)
    flank_3_seq = SL_RNA_SEQUENCE[ASO_TARGET_END:]
    flank_3 = build_single_strand_rna(flank_3_seq, chain_id="A")
    # Ajustar numeracao e posicao z
    z_offset_3 = ASO_TARGET_END * config.rise_per_bp
    for atom in flank_3:
        atom.residue_num += ASO_TARGET_END
        atom.z += z_offset_3
    complex_atoms.extend(flank_3)

    return {
        "sl_rna_free": sl_rna_free,
        "duplex": duplex,
        "complex": complex_atoms,
    }
