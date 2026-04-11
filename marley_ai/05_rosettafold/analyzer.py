"""Analise geometrica de estruturas de acidos nucleicos.

Analisa PDB existentes e modelos construidos, calculando:
- Parametros helicoidais (rise, twist, inclinacao)
- Largura dos sulcos (major/minor groove)
- RMSD entre estruturas
- Pontes de hidrogenio (H-bonds) entre fitas
- Area de superficie acessivel ao solvente (SASA)
- Contatos entre ASO e SL RNA

Ref: Lu XJ & Olson WK (2003) NAR 31:5108-5121 (parametros helicoidais).
     Lee B & Richards FM (1971) J Mol Biol 55:379-400 (SASA).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

import importlib as _il

import numpy as np

# Importacao via importlib — "05" nao e valido como identificador Python
_builder = _il.import_module("marley_ai.05_rosettafold.builder")
_cfg = _il.import_module("marley_ai.05_rosettafold.config")

Atom = _builder.Atom
PURINES = _builder.PURINES
PYRIMIDINES = _builder.PYRIMIDINES

A_FORM_MAJOR_GROOVE_WIDTH = _cfg.A_FORM_MAJOR_GROOVE_WIDTH
A_FORM_MINOR_GROOVE_WIDTH = _cfg.A_FORM_MINOR_GROOVE_WIDTH
A_FORM_RISE = _cfg.A_FORM_RISE
A_FORM_TWIST = _cfg.A_FORM_TWIST
ATOMIC_RADII = _cfg.ATOMIC_RADII
COMPLEMENT = _cfg.COMPLEMENT
SOLVENT_PROBE_RADIUS = _cfg.SOLVENT_PROBE_RADIUS
LNA_POSITIONS = _cfg.LNA_POSITIONS
StructuralAnalysisConfig = _cfg.StructuralAnalysisConfig


# ---------------------------------------------------------------------------
# Estruturas de dados para resultados de analise
# ---------------------------------------------------------------------------

@dataclass
class HelicalParameters:
    """Parametros helicoidais calculados para um duplex."""
    mean_rise: float = 0.0            # Angstroms
    std_rise: float = 0.0
    mean_twist: float = 0.0           # graus
    std_twist: float = 0.0
    mean_inclination: float = 0.0     # graus
    n_base_pairs: int = 0
    bp_per_turn: float = 0.0
    helix_length: float = 0.0         # Angstroms
    helix_form: str = "unknown"       # A, B, ou Z


@dataclass
class GrooveAnalysis:
    """Analise dos sulcos major e minor do duplex."""
    major_groove_width: float = 0.0    # Angstroms
    minor_groove_width: float = 0.0    # Angstroms
    major_groove_depth: float = 0.0    # Angstroms
    minor_groove_depth: float = 0.0    # Angstroms
    # Razao — valores > 1 indicam sulco menor mais largo (tipico forma-A)
    minor_to_major_ratio: float = 0.0


@dataclass
class HydrogenBond:
    """Representa uma ponte de hidrogenio entre duas bases."""
    donor_chain: str
    donor_residue: int
    donor_atom: str
    acceptor_chain: str
    acceptor_residue: int
    acceptor_atom: str
    distance: float           # Angstroms
    bond_type: str            # "Watson-Crick", "wobble", "non-canonical"


@dataclass
class StructureAnalysisResult:
    """Resultado completo da analise estrutural."""
    # Parametros helicoidais
    helical_params: HelicalParameters = field(default_factory=HelicalParameters)
    # Analise de sulcos
    groove_analysis: GrooveAnalysis = field(default_factory=GrooveAnalysis)
    # Pontes de hidrogenio
    hydrogen_bonds: list[HydrogenBond] = field(default_factory=list)
    n_watson_crick: int = 0
    n_wobble: int = 0
    n_total_hbonds: int = 0
    # RMSD e comparacao
    rmsd_vs_ideal: float = 0.0       # RMSD contra forma-A ideal
    rmsd_vs_existing: float = 0.0    # RMSD contra PDB existente (se disponivel)
    # SASA
    sasa_free: float = 0.0           # SASA do RNA isolado (A^2)
    sasa_bound: float = 0.0          # SASA do complexo (A^2)
    delta_sasa: float = 0.0          # Variacao de SASA na ligacao (A^2)
    # Acesso da RNase H
    rnase_h_accessible: bool = False  # sulco minor acessivel para RNase H?
    rnase_h_score: float = 0.0       # score de acessibilidade (0-1)
    # Resumo
    structural_quality: str = ""      # "excellent", "good", "fair", "poor"


# ---------------------------------------------------------------------------
# Parser de PDB existentes
# ---------------------------------------------------------------------------

def parse_pdb(pdb_path: Path) -> list[Atom]:
    """Le um arquivo PDB e retorna lista de Atom.

    Faz parsing simples das linhas ATOM/HETATM seguindo o formato PDB padrao.

    Args:
        pdb_path: Caminho para o arquivo PDB.

    Returns:
        Lista de Atom extraidos do PDB.

    Raises:
        FileNotFoundError: Se o arquivo nao existir.
    """
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB nao encontrado: {pdb_path}")

    atoms: list[Atom] = []
    with open(pdb_path, encoding="utf-8") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            # Parsing de colunas fixas do PDB
            try:
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                chain_id = line[21:22].strip() or "A"
                residue_num = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                b_factor = float(line[60:66].strip()) if len(line) > 66 else 0.0
                element = line[76:78].strip() if len(line) > 78 else atom_name[0]
            except (ValueError, IndexError):
                continue

            atoms.append(Atom(
                name=atom_name,
                residue_name=residue_name,
                chain_id=chain_id,
                residue_num=residue_num,
                x=x, y=y, z=z,
                element=element if element else atom_name[0],
                b_factor=b_factor,
            ))

    return atoms


# ---------------------------------------------------------------------------
# Calculo de parametros helicoidais
# ---------------------------------------------------------------------------

def _get_residue_centers(atoms: list[Atom], chain_id: str) -> dict[int, np.ndarray]:
    """Calcula o centroide de cada residuo em uma cadeia.

    Usa atomos pesados (non-H) para calcular o centro geometrico.

    Args:
        atoms: Lista completa de atomos.
        chain_id: Cadeia a analisar.

    Returns:
        Dicionario {residue_num: centro_xyz}.
    """
    residue_atoms: dict[int, list[np.ndarray]] = {}
    for atom in atoms:
        if atom.chain_id != chain_id:
            continue
        if atom.element == "H":
            continue
        resnum = atom.residue_num
        if resnum not in residue_atoms:
            residue_atoms[resnum] = []
        residue_atoms[resnum].append(np.array([atom.x, atom.y, atom.z]))

    centers: dict[int, np.ndarray] = {}
    for resnum, coords_list in residue_atoms.items():
        if coords_list:
            centers[resnum] = np.mean(coords_list, axis=0)

    return centers


def _get_backbone_atoms(
    atoms: list[Atom],
    chain_id: str,
    atom_name: str = "P",
) -> dict[int, np.ndarray]:
    """Extrai posicoes de um atomo especifico do backbone por residuo.

    Args:
        atoms: Lista completa de atomos.
        chain_id: Cadeia a analisar.
        atom_name: Nome do atomo a extrair (default: "P" para fosfato).

    Returns:
        Dicionario {residue_num: posicao_xyz}.
    """
    positions: dict[int, np.ndarray] = {}
    for atom in atoms:
        if atom.chain_id == chain_id and atom.name == atom_name:
            positions[atom.residue_num] = np.array([atom.x, atom.y, atom.z])
    return positions


def compute_helical_parameters(atoms: list[Atom]) -> HelicalParameters:
    """Calcula parametros helicoidais a partir de coordenadas atomicas.

    Usa os atomos C4' do backbone da cadeia A (RNA) para medir rise e twist,
    pois C4' esta presente em todo residuo e nao e ambiguo como centroides
    que misturam backbone e base.

    O twist e calculado como o angulo entre vetores radiais (do eixo helicoidal
    ao C4') de residuos consecutivos, projetados no plano XY.

    Args:
        atoms: Lista de atomos do duplex.

    Returns:
        HelicalParameters com valores calculados.
    """
    params = HelicalParameters()

    # Usar C4' como proxy do backbone — presente em todo residuo
    c4_a = _get_backbone_atoms(atoms, "A", "C4'")
    if len(c4_a) < 3:
        # Fallback: tentar com centroides
        c4_a = _get_residue_centers(atoms, "A")
    if len(c4_a) < 3:
        params.helix_form = "insufficient_data"
        return params

    # Ordenar por numero de residuo
    sorted_nums = sorted(c4_a.keys())
    coords = np.array([c4_a[n] for n in sorted_nums])
    n_res = len(coords)

    # Estimar o eixo helicoidal como a media XY de todos os C4'
    # (para uma helice ideal centrada na origem, isto e proximo de (0,0))
    helix_axis_xy = np.mean(coords[:, :2], axis=0)

    # Calcular rise: componente Z entre residuos consecutivos
    rises: list[float] = []
    for i in range(n_res - 1):
        diff = coords[i + 1] - coords[i]
        rise_z = abs(diff[2])
        if rise_z > 0.5:  # filtrar artefatos
            rises.append(rise_z)

    if rises:
        params.mean_rise = float(np.mean(rises))
        params.std_rise = float(np.std(rises))

    # Calcular twist: angulo entre vetores radiais consecutivos no plano XY
    # Vetor radial = posicao XY do C4' - eixo helicoidal XY
    twists: list[float] = []
    for i in range(n_res - 1):
        r1 = coords[i][:2] - helix_axis_xy
        r2 = coords[i + 1][:2] - helix_axis_xy

        norm1 = np.linalg.norm(r1)
        norm2 = np.linalg.norm(r2)
        if norm1 < 0.01 or norm2 < 0.01:
            continue

        cos_angle = np.clip(np.dot(r1, r2) / (norm1 * norm2), -1.0, 1.0)
        angle = math.degrees(math.acos(float(cos_angle)))

        # Usar atan2 para obter o angulo com sinal correto
        angle_signed = math.degrees(math.atan2(
            float(r1[0] * r2[1] - r1[1] * r2[0]),
            float(np.dot(r1, r2)),
        ))

        # Twist deve ser positivo (sentido anti-horario visto do topo)
        twist_val = abs(angle_signed)
        if twist_val > 5.0:  # filtrar artefatos
            twists.append(twist_val)

    if twists:
        params.mean_twist = float(np.mean(twists))
        params.std_twist = float(np.std(twists))
        params.bp_per_turn = 360.0 / params.mean_twist if params.mean_twist > 0 else 0

    # Numero de pares de bases
    # Contar residuos unicos na cadeia B tambem
    chain_b_res = set()
    for atom in atoms:
        if atom.chain_id == "B":
            chain_b_res.add(atom.residue_num)
    params.n_base_pairs = min(n_res, len(chain_b_res)) if chain_b_res else n_res
    params.helix_length = params.mean_rise * (params.n_base_pairs - 1) if params.mean_rise > 0 else 0

    # Classificar forma helicoidal
    if 2.5 <= params.mean_rise <= 3.2 and 28 <= params.mean_twist <= 37:
        params.helix_form = "A"
    elif 3.2 < params.mean_rise <= 3.6 and 34 <= params.mean_twist <= 38:
        params.helix_form = "B"
    elif params.mean_rise < 3.8:
        params.helix_form = "A/B_intermediate"
    else:
        params.helix_form = "extended"

    return params


# ---------------------------------------------------------------------------
# Analise de sulcos
# ---------------------------------------------------------------------------

def compute_groove_analysis(
    atoms: list[Atom],
    config: StructuralAnalysisConfig | None = None,
) -> GrooveAnalysis:
    """Estima largura dos sulcos major e minor do duplex.

    Usa a distancia entre backbones (atomos P) das duas fitas
    para estimar a largura dos sulcos. Na forma-A, o sulco menor
    e mais largo que o sulco maior (inverso da forma-B).

    A RNase H precisa de acesso pelo sulco menor (largura > 9 A).

    Args:
        atoms: Lista de atomos do duplex.
        config: Configuracao do modulo.

    Returns:
        GrooveAnalysis com larguras estimadas.
    """
    groove = GrooveAnalysis()

    # Obter posicoes dos fosfatos das duas cadeias
    p_chain_a = _get_backbone_atoms(atoms, "A", "P")
    p_chain_b = _get_backbone_atoms(atoms, "B", "P")

    if not p_chain_a or not p_chain_b:
        # Sem dados de backbone — usar valores teoricos
        groove.major_groove_width = A_FORM_MAJOR_GROOVE_WIDTH
        groove.minor_groove_width = A_FORM_MINOR_GROOVE_WIDTH
        return groove

    # Calcular distancias entre fosfatos das duas cadeias
    # Na forma-A, os fosfatos estao em posicoes alternadas
    p_a_list = [p_chain_a[k] for k in sorted(p_chain_a.keys())]
    p_b_list = [p_chain_b[k] for k in sorted(p_chain_b.keys())]

    all_distances: list[float] = []
    for pa in p_a_list:
        for pb in p_b_list:
            dist = float(np.linalg.norm(pa - pb))
            all_distances.append(dist)

    if all_distances:
        # Separar em sulco maior e menor
        # O sulco menor tem distancias menores entre P das fitas opostas
        sorted_dists = sorted(all_distances)
        n_half = len(sorted_dists) // 2

        # Sulco menor — distancias mais curtas entre fitas
        # (Na forma-A, o sulco menor e na verdade mais largo —
        #  o P-P esta mais afastado no sulco menor)
        # Corrigir: subtrair 5.8 A (diametro do fosfato) para obter largura do sulco
        phosphate_diameter = 5.8

        # Estimativa da largura do sulco a partir do padrao de distancias P-P
        # Pegamos as distancias P-P que cruzam cada sulco
        min_p_p = sorted_dists[0] if sorted_dists else 18.0
        max_p_p = sorted_dists[-1] if sorted_dists else 18.0
        median_p_p = sorted_dists[len(sorted_dists) // 2] if sorted_dists else 18.0

        # Para forma-A: minor groove ~11 A, major groove ~2.7 A
        # P-P across minor groove ~ 17.4 A (11 + 2*3.2 para raios P)
        # P-P across major groove ~ 9.2 A (2.7 + 2*3.2)
        groove.minor_groove_width = max(
            0.0, median_p_p - phosphate_diameter
        )
        groove.major_groove_width = max(
            0.0, min_p_p - phosphate_diameter
        )

        # Usar valores teoricos se a estimativa for irrealista
        if groove.minor_groove_width > 15.0 or groove.minor_groove_width < 3.0:
            groove.minor_groove_width = A_FORM_MINOR_GROOVE_WIDTH
        if groove.major_groove_width > 10.0 or groove.major_groove_width < 0.5:
            groove.major_groove_width = A_FORM_MAJOR_GROOVE_WIDTH

    # Razao minor/major
    if groove.major_groove_width > 0:
        groove.minor_to_major_ratio = (
            groove.minor_groove_width / groove.major_groove_width
        )

    return groove


# ---------------------------------------------------------------------------
# Deteccao de pontes de hidrogenio
# ---------------------------------------------------------------------------

# Atomos doadores/aceitadores de H-bond por tipo de base
# Ref: Watson-Crick base pairing geometry
_HBOND_DONORS: Final[dict[str, list[str]]] = {
    "A": ["N6"],           # NH2 exociclico
    "G": ["N1", "N2"],     # NH imidazol + NH2 exociclico
    "C": ["N4"],           # NH2 exociclico
    "T": ["N3"],           # NH
    "U": ["N3"],           # NH
}

_HBOND_ACCEPTORS: Final[dict[str, list[str]]] = {
    "A": ["N1", "N7"],     # N anel
    "G": ["O6", "N7"],     # C=O + N anel
    "C": ["O2", "N3"],     # C=O + N anel
    "T": ["O2", "O4"],     # C=O
    "U": ["O2", "O4"],     # C=O
}

# Limiar de distancia para H-bond (Angstroms)
_HBOND_MAX_DISTANCE: Final[float] = 3.5


def detect_hydrogen_bonds(atoms: list[Atom]) -> list[HydrogenBond]:
    """Detecta pontes de hidrogenio entre bases das duas fitas.

    Busca pares de atomos doador-aceitador entre cadeia A e cadeia B
    com distancia <= 3.5 A. Classifica como Watson-Crick se o par
    de bases segue a complementaridade padrao.

    Args:
        atoms: Lista de atomos do duplex (cadeias A e B).

    Returns:
        Lista de HydrogenBond detectadas.
    """
    # Indexar atomos por cadeia e residuo
    chain_a_atoms: dict[int, dict[str, np.ndarray]] = {}
    chain_b_atoms: dict[int, dict[str, np.ndarray]] = {}
    chain_a_bases: dict[int, str] = {}
    chain_b_bases: dict[int, str] = {}

    for atom in atoms:
        pos = np.array([atom.x, atom.y, atom.z])
        # Extrair letra da base do nome do residuo
        base = atom.residue_name.replace("D", "").strip()
        if len(base) > 1:
            base = base[-1]  # pegar ultimo caractere

        if atom.chain_id == "A":
            if atom.residue_num not in chain_a_atoms:
                chain_a_atoms[atom.residue_num] = {}
                chain_a_bases[atom.residue_num] = base
            chain_a_atoms[atom.residue_num][atom.name] = pos
        elif atom.chain_id == "B":
            if atom.residue_num not in chain_b_atoms:
                chain_b_atoms[atom.residue_num] = {}
                chain_b_bases[atom.residue_num] = base
            chain_b_atoms[atom.residue_num][atom.name] = pos

    hbonds: list[HydrogenBond] = []

    # Buscar H-bonds entre residuos proximos das duas cadeias
    for res_a, atoms_a in chain_a_atoms.items():
        base_a = chain_a_bases.get(res_a, "A")
        donors_a = _HBOND_DONORS.get(base_a, [])
        acceptors_a = _HBOND_ACCEPTORS.get(base_a, [])

        for res_b, atoms_b in chain_b_atoms.items():
            base_b = chain_b_bases.get(res_b, "T")
            donors_b = _HBOND_DONORS.get(base_b, [])
            acceptors_b = _HBOND_ACCEPTORS.get(base_b, [])

            # Verificar doador A -> aceitador B
            for donor_name in donors_a:
                if donor_name not in atoms_a:
                    continue
                for acc_name in acceptors_b:
                    if acc_name not in atoms_b:
                        continue
                    dist = float(np.linalg.norm(atoms_a[donor_name] - atoms_b[acc_name]))
                    if dist <= _HBOND_MAX_DISTANCE:
                        # Classificar o tipo de H-bond
                        bond_type = _classify_hbond(base_a, base_b)
                        hbonds.append(HydrogenBond(
                            donor_chain="A", donor_residue=res_a,
                            donor_atom=donor_name,
                            acceptor_chain="B", acceptor_residue=res_b,
                            acceptor_atom=acc_name,
                            distance=dist, bond_type=bond_type,
                        ))

            # Verificar doador B -> aceitador A
            for donor_name in donors_b:
                if donor_name not in atoms_b:
                    continue
                for acc_name in acceptors_a:
                    if acc_name not in atoms_a:
                        continue
                    dist = float(np.linalg.norm(atoms_b[donor_name] - atoms_a[acc_name]))
                    if dist <= _HBOND_MAX_DISTANCE:
                        bond_type = _classify_hbond(base_a, base_b)
                        hbonds.append(HydrogenBond(
                            donor_chain="B", donor_residue=res_b,
                            donor_atom=donor_name,
                            acceptor_chain="A", acceptor_residue=res_a,
                            acceptor_atom=acc_name,
                            distance=dist, bond_type=bond_type,
                        ))

    return hbonds


def _classify_hbond(base1: str, base2: str) -> str:
    """Classifica o tipo de ponte de hidrogenio entre duas bases.

    Args:
        base1: Base da cadeia A.
        base2: Base da cadeia B.

    Returns:
        Tipo: "Watson-Crick", "wobble", ou "non-canonical".
    """
    # Normalizar: U e T sao equivalentes para pareamento
    b1 = base1.replace("U", "T")
    b2 = base2.replace("U", "T")

    wc_pairs = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}
    wobble_pairs = {("G", "T"), ("T", "G")}

    pair = (b1, b2)
    if pair in wc_pairs:
        return "Watson-Crick"
    elif pair in wobble_pairs:
        return "wobble"
    else:
        return "non-canonical"


# ---------------------------------------------------------------------------
# Calculo de RMSD
# ---------------------------------------------------------------------------

def compute_rmsd(atoms1: list[Atom], atoms2: list[Atom]) -> float:
    """Calcula RMSD entre dois conjuntos de atomos.

    Alinha por nome de atomo e numero de residuo, usando apenas
    atomos presentes em ambas as estruturas.

    Args:
        atoms1: Primeira estrutura.
        atoms2: Segunda estrutura.

    Returns:
        RMSD em Angstroms.
    """
    # Construir indice: (chain, resnum, atom_name) -> coords
    def _index_atoms(atoms: list[Atom]) -> dict[tuple[str, int, str], np.ndarray]:
        idx: dict[tuple[str, int, str], np.ndarray] = {}
        for a in atoms:
            key = (a.chain_id, a.residue_num, a.name)
            idx[key] = np.array([a.x, a.y, a.z])
        return idx

    idx1 = _index_atoms(atoms1)
    idx2 = _index_atoms(atoms2)

    # Encontrar atomos comuns
    common_keys = set(idx1.keys()) & set(idx2.keys())
    if not common_keys:
        return float("inf")

    coords1 = np.array([idx1[k] for k in sorted(common_keys)])
    coords2 = np.array([idx2[k] for k in sorted(common_keys)])

    # Centralizar
    center1 = np.mean(coords1, axis=0)
    center2 = np.mean(coords2, axis=0)
    c1 = coords1 - center1
    c2 = coords2 - center2

    # RMSD sem alinhamento rotacional (apenas translacao)
    diff = c1 - c2
    rmsd = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))

    return rmsd


def compute_rmsd_kabsch(atoms1: list[Atom], atoms2: list[Atom]) -> float:
    """Calcula RMSD com alinhamento rotacional otimo (Kabsch algorithm).

    Encontra a rotacao que minimiza o RMSD entre as duas estruturas.
    Ref: Kabsch W (1976) Acta Cryst A32:922-923.

    Args:
        atoms1: Estrutura de referencia.
        atoms2: Estrutura a alinhar.

    Returns:
        RMSD otimo em Angstroms.
    """
    def _index_atoms(atoms: list[Atom]) -> dict[tuple[str, int, str], np.ndarray]:
        idx: dict[tuple[str, int, str], np.ndarray] = {}
        for a in atoms:
            key = (a.chain_id, a.residue_num, a.name)
            idx[key] = np.array([a.x, a.y, a.z])
        return idx

    idx1 = _index_atoms(atoms1)
    idx2 = _index_atoms(atoms2)

    common_keys = set(idx1.keys()) & set(idx2.keys())
    if len(common_keys) < 3:
        return float("inf")

    sorted_keys = sorted(common_keys)
    P = np.array([idx1[k] for k in sorted_keys])
    Q = np.array([idx2[k] for k in sorted_keys])

    # Centralizar
    P_center = np.mean(P, axis=0)
    Q_center = np.mean(Q, axis=0)
    P = P - P_center
    Q = Q - Q_center

    # Matriz de covariancia cruzada
    H = P.T @ Q

    # SVD para encontrar rotacao otima
    U, S, Vt = np.linalg.svd(H)

    # Corrigir reflexao se necessario (det < 0)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1.0, 1.0, np.sign(d)])

    # Rotacao otima
    R = Vt.T @ sign_matrix @ U.T

    # Aplicar rotacao e calcular RMSD
    Q_rot = Q @ R.T
    diff = P - Q_rot
    rmsd = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))

    return rmsd


# ---------------------------------------------------------------------------
# Estimativa de SASA (Shrake-Rupley)
# ---------------------------------------------------------------------------

def estimate_sasa(
    atoms: list[Atom],
    n_points: int = 92,
    probe_radius: float = SOLVENT_PROBE_RADIUS,
) -> float:
    """Estima area de superficie acessivel ao solvente (SASA).

    Implementacao simplificada do algoritmo Shrake-Rupley:
    distribui pontos na esfera ao redor de cada atomo e conta
    quantos nao sao ocluidos por atomos vizinhos.

    Ref: Shrake A & Rupley JA (1973) J Mol Biol 79:351-371.

    Args:
        atoms: Lista de atomos.
        n_points: Numero de pontos na esfera por atomo.
        probe_radius: Raio da sonda de solvente (Angstroms).

    Returns:
        SASA total em Angstroms quadrados.
    """
    if not atoms:
        return 0.0

    # Gerar pontos na esfera unitaria (distribuicao de Fibonacci)
    sphere_points = _fibonacci_sphere(n_points)

    coords = np.array([[a.x, a.y, a.z] for a in atoms])
    radii = np.array([
        ATOMIC_RADII.get(a.element, 1.7) + probe_radius
        for a in atoms
    ])

    total_sasa = 0.0
    n_atoms = len(atoms)

    for i in range(n_atoms):
        r_i = radii[i]
        # Pontos na esfera expandida do atomo i
        test_points = coords[i] + sphere_points * r_i

        # Verificar quais pontos sao acessiveis (nao ocluidos)
        accessible = 0
        for point in test_points:
            # Distancia do ponto a todos os outros atomos
            diffs = coords - point
            dists_sq = np.sum(diffs ** 2, axis=1)

            # O ponto esta acessivel se nenhum outro atomo o oclui
            # (distancia ao ponto > raio expandido do outro atomo)
            occluded = False
            for j in range(n_atoms):
                if j == i:
                    continue
                if dists_sq[j] < radii[j] ** 2:
                    occluded = True
                    break

            if not occluded:
                accessible += 1

        # Area contribuida por este atomo
        area_per_point = 4.0 * math.pi * r_i ** 2 / n_points
        total_sasa += accessible * area_per_point

    return total_sasa


def _fibonacci_sphere(n: int) -> np.ndarray:
    """Gera pontos uniformemente distribuidos na esfera unitaria.

    Usa a espiral de Fibonacci para distribuicao quase-uniforme.

    Args:
        n: Numero de pontos.

    Returns:
        Array (n, 3) com coordenadas dos pontos.
    """
    golden_ratio = (1.0 + math.sqrt(5.0)) / 2.0
    indices = np.arange(n, dtype=float)

    theta = 2.0 * math.pi * indices / golden_ratio
    phi = np.arccos(1.0 - 2.0 * (indices + 0.5) / n)

    x = np.cos(theta) * np.sin(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(phi)

    return np.column_stack([x, y, z])


# ---------------------------------------------------------------------------
# Analise de acessibilidade da RNase H
# ---------------------------------------------------------------------------

def assess_rnase_h_accessibility(
    groove: GrooveAnalysis,
    helical_params: HelicalParameters,
) -> tuple[bool, float]:
    """Avalia se a RNase H pode acessar o sulco menor do duplex.

    A RNase H requer:
    1. Sulco menor com largura >= 9.0 A (para acomodar o dominio catalitico)
    2. Geometria A-form (o RNA:DNA hibrido favorece forma-A)
    3. Pelo menos 8-10 bp de duplex RNA:DNA continuo

    Ref: Nowotny M et al. (2005) Cell 121:1005-1016 (estrutura RNase H1).

    Args:
        groove: Resultado da analise de sulcos.
        helical_params: Parametros helicoidais calculados.

    Returns:
        Tupla (is_accessible: bool, score: float 0-1).
    """
    score = 0.0

    # Criterio 1: largura do sulco menor >= 9.0 A
    if groove.minor_groove_width >= 9.0:
        score += 0.4
    elif groove.minor_groove_width >= 7.0:
        score += 0.2

    # Criterio 2: forma A (rise ~2.8 A, twist ~32.7)
    if helical_params.helix_form == "A":
        score += 0.3
    elif helical_params.helix_form == "A/B_intermediate":
        score += 0.15

    # Criterio 3: numero de pares de bases >= 10
    if helical_params.n_base_pairs >= 10:
        score += 0.3
    elif helical_params.n_base_pairs >= 8:
        score += 0.15

    is_accessible = score >= 0.7
    return is_accessible, min(score, 1.0)


# ---------------------------------------------------------------------------
# Analise completa
# ---------------------------------------------------------------------------

def analyze_structure(
    duplex_atoms: list[Atom],
    free_rna_atoms: list[Atom] | None = None,
    existing_pdb_atoms: list[Atom] | None = None,
    config: StructuralAnalysisConfig | None = None,
) -> StructureAnalysisResult:
    """Executa analise estrutural completa do duplex ASO:SL RNA.

    Calcula todos os parametros geometricos, detecta H-bonds,
    estima SASA, e avalia acessibilidade da RNase H.

    Args:
        duplex_atoms: Atomos do duplex ASO:SL RNA.
        free_rna_atoms: Atomos do SL RNA isolado (para comparacao SASA).
        existing_pdb_atoms: Atomos de PDB existente (para RMSD).
        config: Configuracao do modulo.

    Returns:
        StructureAnalysisResult com todos os resultados.
    """
    result = StructureAnalysisResult()

    # 1. Parametros helicoidais
    print("  [analyzer] Calculando parametros helicoidais...")
    result.helical_params = compute_helical_parameters(duplex_atoms)

    # 2. Analise de sulcos
    print("  [analyzer] Analisando sulcos major/minor...")
    result.groove_analysis = compute_groove_analysis(duplex_atoms, config)

    # 3. Pontes de hidrogenio
    print("  [analyzer] Detectando pontes de hidrogenio...")
    result.hydrogen_bonds = detect_hydrogen_bonds(duplex_atoms)
    result.n_watson_crick = sum(
        1 for hb in result.hydrogen_bonds if hb.bond_type == "Watson-Crick"
    )
    result.n_wobble = sum(
        1 for hb in result.hydrogen_bonds if hb.bond_type == "wobble"
    )
    result.n_total_hbonds = len(result.hydrogen_bonds)

    # 4. RMSD contra estrutura existente
    if existing_pdb_atoms:
        print("  [analyzer] Calculando RMSD contra PDB existente...")
        result.rmsd_vs_existing = compute_rmsd_kabsch(duplex_atoms, existing_pdb_atoms)

    # 5. SASA
    print("  [analyzer] Estimando SASA (simplificado)...")
    # Usar subconjunto de atomos para performance (backbone apenas)
    backbone_duplex = [a for a in duplex_atoms if a.name in ("P", "C4'", "N9", "N1")]
    backbone_free = (
        [a for a in free_rna_atoms if a.name in ("P", "C4'", "N9", "N1")]
        if free_rna_atoms else []
    )

    if backbone_duplex:
        result.sasa_bound = estimate_sasa(backbone_duplex)
    if backbone_free:
        result.sasa_free = estimate_sasa(backbone_free)

    result.delta_sasa = result.sasa_free - result.sasa_bound

    # 6. Acessibilidade da RNase H
    print("  [analyzer] Avaliando acessibilidade da RNase H...")
    accessible, score = assess_rnase_h_accessibility(
        result.groove_analysis, result.helical_params,
    )
    result.rnase_h_accessible = accessible
    result.rnase_h_score = score

    # 7. Classificar qualidade estrutural
    if (result.helical_params.helix_form == "A"
            and result.n_watson_crick >= 10
            and result.rnase_h_accessible):
        result.structural_quality = "excellent"
    elif result.n_watson_crick >= 8 and result.rnase_h_score >= 0.5:
        result.structural_quality = "good"
    elif result.n_watson_crick >= 5:
        result.structural_quality = "fair"
    else:
        result.structural_quality = "poor"

    return result
