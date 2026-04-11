"""Extracao de nuvens de pontos 3D a partir de estruturas PDB.

Estrategia:
- Parseia arquivo PDB e extrai coordenadas de atomos do backbone
  (P, C3', C4', O3', O5') que definem a geometria da cadeia principal.
- Se o PDB nao existe ou e incompleto, gera backbone sintetico usando
  geometria B-DNA idealizada (rise 3.38 A, twist 36 graus, raio 10 A).
- Produz tres nuvens de pontos: ASO livre, SL RNA alvo, duplex ASO:SL RNA.

Referencia geometrica:
    Watson & Crick (1953); Saenger W (1984) "Principles of Nucleic Acid Structure"
    Rise por par de bases: 3.38 A (B-DNA canonico)
    Twist por par de bases: 36 graus (10 bp por volta completa)
    Raio do backbone sugar-fosfato: ~10 A
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from core.logger import get_logger

logger = get_logger("math_3_tda.point_cloud")

# ---------------------------------------------------------------------------
# Constantes geometricas de B-DNA idealizado
# ---------------------------------------------------------------------------

# Subida por par de bases ao longo do eixo helicoidal (Angstroms)
RISE_PER_BP: float = 3.38

# Rotacao por par de bases (graus) — 10 bp completam uma volta
TWIST_PER_BP_DEG: float = 36.0

# Raio do backbone sugar-fosfato (distancia do eixo helicoidal, Angstroms)
BACKBONE_RADIUS: float = 10.0

# Atomos do backbone que definem a geometria da cadeia principal
BACKBONE_ATOMS: frozenset[str] = frozenset({"P", "C3'", "C4'", "O3'", "O5'"})

# Deslocamento angular entre as duas fitas do duplex (graus)
# Na B-DNA, os sulcos menor e maior dividem os 360 graus de forma assimetrica.
# O sulco menor corresponde a ~150 graus entre as fitas.
STRAND_OFFSET_DEG: float = 150.0


# ---------------------------------------------------------------------------
# Parser de PDB
# ---------------------------------------------------------------------------


def parse_pdb_backbone(pdb_path: Path) -> np.ndarray:
    """Extrai coordenadas (x, y, z) de atomos do backbone de um arquivo PDB.

    Filtra apenas atomos P, C3', C4', O3', O5' de registros ATOM/HETATM.
    O formato PDB (v3.3) define colunas fixas para nome do atomo e coordenadas.

    Args:
        pdb_path: Caminho para o arquivo PDB.

    Returns:
        Array Nx3 de coordenadas em Angstroms. Vazio se arquivo nao existe.
    """
    if not pdb_path.is_file():
        logger.warning("Arquivo PDB nao encontrado: %s", pdb_path)
        return np.empty((0, 3), dtype=np.float64)

    coords: list[list[float]] = []

    with open(pdb_path, encoding="utf-8") as fh:
        for line in fh:
            record = line[:6].strip()
            if record not in ("ATOM", "HETATM"):
                continue

            # Colunas PDB v3.3: nome do atomo em 13-16 (0-indexed: 12-15)
            atom_name = line[12:16].strip()

            if atom_name not in BACKBONE_ATOMS:
                continue

            # Coordenadas: x(31-38), y(39-46), z(47-54)
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
            except (ValueError, IndexError):
                # Linha malformada — pula sem interromper
                continue

    if not coords:
        logger.warning("Nenhum atomo de backbone encontrado em: %s", pdb_path)
        return np.empty((0, 3), dtype=np.float64)

    result = np.array(coords, dtype=np.float64)
    logger.info(
        "PDB %s: %d atomos de backbone extraidos",
        pdb_path.name, len(result),
    )
    return result


# ---------------------------------------------------------------------------
# Gerador de backbone sintetico (B-DNA idealizado)
# ---------------------------------------------------------------------------


def generate_synthetic_backbone(
    sequence: str,
    strand_offset_deg: float = 0.0,
    z_offset: float = 0.0,
) -> np.ndarray:
    """Gera coordenadas de backbone sintetico usando geometria B-DNA idealizada.

    Cada nucleotideo contribui com 5 atomos de backbone (P, C3', C4', O3', O5')
    posicionados ao longo de uma helice com parametros canonicos.

    A posicao de cada atomo e calculada por:
        x = R * cos(theta + delta_i)
        y = R * sin(theta + delta_i)
        z = n * rise + z_offset_i

    onde delta_i e z_offset_i sao deslocamentos locais para cada tipo de atomo,
    baseados na geometria media de B-DNA.

    Args:
        sequence: Sequencia de nucleotideos (ex: "ACAGAAACTG...").
        strand_offset_deg: Deslocamento angular para a segunda fita (graus).
        z_offset: Deslocamento vertical para posicionar a fita (Angstroms).

    Returns:
        Array Nx3 de coordenadas sinteticas (N = 5 * len(sequence)).
    """
    n_bases = len(sequence)
    twist_rad = math.radians(TWIST_PER_BP_DEG)
    strand_offset_rad = math.radians(strand_offset_deg)

    # Deslocamentos locais de cada atomo de backbone relativo ao fosforo
    # Valores aproximados de cristalografia de B-DNA (Saenger, 1984)
    # (delta_angulo em radianos, delta_z em Angstroms, delta_raio em Angstroms)
    atom_offsets: list[tuple[str, float, float, float]] = [
        ("P",   0.0,   0.0,   0.0),     # Fosforo — referencia
        ("O5'", 0.05,  0.3,  -0.5),     # O5' — ligeiramente acima do P
        ("C4'", 0.15,  0.8,  -1.5),     # C4' — parte do acucar
        ("C3'", 0.20,  1.2,  -1.8),     # C3' — acucar desoxirribose
        ("O3'", 0.25,  1.6,  -1.2),     # O3' — ligacao ao proximo fosfato
    ]

    coords: list[list[float]] = []

    for i in range(n_bases):
        # Angulo base deste nucleotideo na helice
        theta_base = i * twist_rad + strand_offset_rad
        z_base = i * RISE_PER_BP + z_offset

        for _atom_name, d_angle, d_z, d_radius in atom_offsets:
            theta = theta_base + d_angle
            r = BACKBONE_RADIUS + d_radius
            z = z_base + d_z

            x = r * math.cos(theta)
            y = r * math.sin(theta)
            coords.append([x, y, z])

    result = np.array(coords, dtype=np.float64)
    logger.info(
        "Backbone sintetico gerado: %d bases -> %d atomos",
        n_bases, len(result),
    )
    return result


# ---------------------------------------------------------------------------
# Construtor de nuvens de pontos para as tres estruturas
# ---------------------------------------------------------------------------


def build_point_clouds(
    aso_sequence: str,
    sl_target_sequence: str,
    aso_pdb_path: Path,
    sl_pdb_path: Path,
) -> dict[str, np.ndarray]:
    """Constroi nuvens de pontos para ASO livre, SL RNA alvo e duplex.

    Estrategia:
    1. Tenta extrair coordenadas do PDB real (se disponivel e com backbone)
    2. Se falhar, gera backbone sintetico a partir da sequencia
    3. Para o duplex, combina as duas fitas com deslocamento angular de 150 graus

    Args:
        aso_sequence: Sequencia do ASO (ex: "ACAGAAACTG...").
        sl_target_sequence: Sequencia da regiao alvo do SL RNA.
        aso_pdb_path: Caminho do PDB do ASO.
        sl_pdb_path: Caminho do PDB do SL RNA.

    Returns:
        Dicionario com chaves 'free_aso', 'free_sl_rna', 'duplex',
        cada uma mapeando para um array Nx3 de coordenadas.
    """
    # --- ASO livre ---
    aso_cloud = parse_pdb_backbone(aso_pdb_path)
    aso_source = "pdb"
    if aso_cloud.shape[0] < 5:
        # PDB vazio ou sem backbone suficiente — usa sintetico
        logger.info("Usando backbone sintetico para ASO (PDB insuficiente)")
        aso_cloud = generate_synthetic_backbone(aso_sequence)
        aso_source = "synthetic"

    # --- SL RNA alvo ---
    sl_cloud = parse_pdb_backbone(sl_pdb_path)
    sl_source = "pdb"
    if sl_cloud.shape[0] < 5:
        logger.info("Usando backbone sintetico para SL RNA (PDB insuficiente)")
        sl_cloud = generate_synthetic_backbone(sl_target_sequence)
        sl_source = "synthetic"

    # --- Duplex ASO:SL RNA ---
    # O duplex e modelado como duas fitas helicoidais complementares,
    # deslocadas angularmente por STRAND_OFFSET_DEG (sulco menor da B-DNA).
    # Usa sempre geometria sintetica para garantir alinhamento correto
    # entre as fitas — PDB real pode ter conformacao incompativel.
    sense_strand = generate_synthetic_backbone(
        sl_target_sequence,
        strand_offset_deg=0.0,
        z_offset=0.0,
    )
    antisense_strand = generate_synthetic_backbone(
        aso_sequence,
        strand_offset_deg=STRAND_OFFSET_DEG,
        z_offset=0.0,
    )
    duplex_cloud = np.vstack([sense_strand, antisense_strand])

    logger.info(
        "Nuvens de pontos construidas: ASO=%d (%s), SL=%d (%s), duplex=%d (synthetic)",
        aso_cloud.shape[0], aso_source,
        sl_cloud.shape[0], sl_source,
        duplex_cloud.shape[0],
    )

    return {
        "free_aso": aso_cloud,
        "free_sl_rna": sl_cloud,
        "duplex": duplex_cloud,
        "sources": {
            "free_aso": aso_source,
            "free_sl_rna": sl_source,
            "duplex": "synthetic",
        },
    }
