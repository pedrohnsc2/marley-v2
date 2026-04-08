"""Gerar estrutura 3D do SL RNA de L. infantum para docking molecular.

Este modulo constroi a estrutura tridimensional do Spliced Leader RNA
(39 nt) a partir da sequencia e estrutura secundaria (dot-bracket).
O resultado e um arquivo PDB pronto para conversao a PDBQT e uso
como receptor em docking com AutoDock Vina.

Estrategias de modelagem 3D (em ordem de prioridade):
1. RNAComposer (via API web) -- modelagem por fragmentos
2. Construcao local com Open Babel (SMILES -> 3D -> PDB)
3. Modelo simplificado helicoidal (fallback puro Python)

Todas as abordagens produzem um PDB valido. A qualidade varia:
RNAComposer > Open Babel > modelo helicoidal.

Uso:
    python -m rna_entropy.11_sl_rna_3d
    python -m rna_entropy.11_sl_rna_3d --force
    python -m rna_entropy.11_sl_rna_3d --dry-run
    python -m rna_entropy.11_sl_rna_3d --method helix
"""

from __future__ import annotations

import json
import math
import subprocess
import urllib.request
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Final

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Constantes
# ---------------------------------------------------------------------------

# Sequencia do Spliced Leader RNA de L. infantum (39 nt)
SL_SEQUENCE: Final[str] = "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG"
SL_LENGTH: Final[int] = len(SL_SEQUENCE)

# Conteudo GC
SL_GC_CONTENT: Final[float] = 0.282

# Diretorios de saida
OUTPUT_DIR: Final[Path] = Path("data/sl_rna_3d")
RESULTS_DIR: Final[Path] = Path("results/sl_drug")

# Arquivos de saida
SL_PDB_PATH: Final[Path] = OUTPUT_DIR / "sl_rna.pdb"
SL_PDBQT_PATH: Final[Path] = OUTPUT_DIR / "sl_rna_receptor.pdbqt"
SL_METADATA_PATH: Final[Path] = RESULTS_DIR / "sl_rna_3d_metadata.json"

# Estrutura secundaria do SL RNA (prevista pelo modulo 06/08)
# Para RNA de 39nt com GC=28.2%, estrutura e majoritariamente aberta
# ViennaRNA preve regioes de stem-loop menores
SL_DOT_BRACKET_FALLBACK: Final[str] = "." * SL_LENGTH

# URL do RNAComposer (servidor web para modelagem 3D de RNA)
RNACOMPOSER_URL: Final[str] = "https://rnacomposer.cs.put.poznan.pl"

# Parametros geometricos para modelo helicoidal A-form RNA
# Parametros da helice A do RNA (angstroms e graus)
RISE_PER_BASE: Final[float] = 2.81  # angstroms por base
TWIST_PER_BASE: Final[float] = 32.7  # graus por base
HELIX_RADIUS: Final[float] = 10.0  # raio da helice em angstroms

# Mapeamento de bases RNA para atomos pesados (simplificado)
# Cada base e representada por P, C4', N1/N9 (purina/pirimidina)
RNA_ATOM_TYPES: Final[dict[str, str]] = {
    "A": "N9",  # purina - N9
    "G": "N9",  # purina - N9
    "C": "N1",  # pirimidina - N1
    "U": "N1",  # pirimidina - N1
    "T": "N1",  # pirimidina (DNA) - tratado como U
}

logger = get_logger("sl_rna_3d")


# ---------------------------------------------------------------------------
# Carregamento da estrutura secundaria existente
# ---------------------------------------------------------------------------


def _load_existing_structure() -> str:
    """Carrega a estrutura secundaria do SL RNA dos resultados anteriores.

    Procura o arquivo sl_secondary_structure.json gerado pelo modulo 08.
    Se nao encontrar, retorna o fallback (tudo aberto).

    Returns:
        Notacao dot-bracket da estrutura secundaria.
    """
    structure_json = Path("results/aso/sl_secondary_structure.json")
    if structure_json.exists():
        with open(structure_json, encoding="utf-8") as fh:
            data = json.load(fh)
        dot_bracket = data.get("dot_bracket", SL_DOT_BRACKET_FALLBACK)
        logger.info(
            "Estrutura secundaria carregada de %s: %s",
            structure_json,
            dot_bracket,
        )
        return dot_bracket

    logger.warning(
        "Arquivo de estrutura secundaria nao encontrado em %s. "
        "Usando fallback (todas as posicoes abertas).",
        structure_json,
    )
    return SL_DOT_BRACKET_FALLBACK


# ---------------------------------------------------------------------------
# Metodo 1: RNAComposer (API web)
# ---------------------------------------------------------------------------


def _try_rnacomposer(sequence: str, dot_bracket: str) -> str | None:
    """Tenta obter estrutura 3D via RNAComposer web API.

    RNAComposer usa modelagem por fragmentos baseada em estruturas
    experimentais de RNA do PDB. Produz modelos de alta qualidade
    para RNAs curtos (< 500 nt).

    NOTA: Este metodo requer conexao com a internet e o servidor
    pode estar indisponivel. E usado como primeira opcao por
    produzir a melhor qualidade.

    Args:
        sequence: Sequencia RNA (deve usar U em vez de T).
        dot_bracket: Estrutura secundaria em notacao dot-bracket.

    Returns:
        Conteudo PDB como string, ou None se falhar.
    """
    # Converter T para U (RNAComposer espera RNA)
    rna_seq = sequence.upper().replace("T", "U")

    logger.info(
        "Tentando RNAComposer para %d nt (sequencia: %s...)",
        len(rna_seq),
        rna_seq[:20],
    )

    try:
        # RNAComposer aceita POST com sequencia e estrutura
        # Formato: sequencia\nestrutura
        payload = f"{rna_seq}\n{dot_bracket}"
        data = urllib.parse.urlencode({"sequence": payload}).encode("utf-8")

        req = urllib.request.Request(
            f"{RNACOMPOSER_URL}/api/predict",
            data=data,
            method="POST",
        )
        req.add_header("Content-Type", "application/x-www-form-urlencoded")

        with urllib.request.urlopen(req, timeout=60) as resp:
            if resp.status == 200:
                pdb_content = resp.read().decode("utf-8")
                if "ATOM" in pdb_content:
                    logger.info(
                        "RNAComposer retornou estrutura 3D com sucesso."
                    )
                    return pdb_content
                logger.warning(
                    "RNAComposer retornou resposta sem ATOM records."
                )
                return None
            logger.warning(
                "RNAComposer retornou status %d.", resp.status,
            )
            return None

    except Exception as exc:
        logger.warning(
            "RNAComposer indisponivel: %s. Tentando metodo alternativo.",
            exc,
        )
        return None


# ---------------------------------------------------------------------------
# Metodo 2: Open Babel (construcao local)
# ---------------------------------------------------------------------------


def _try_openbabel_build(sequence: str) -> str | None:
    """Constroi estrutura 3D do RNA usando Open Babel.

    Converte a sequencia RNA para SMILES simplificado de cada
    nucleotideo, gera coordenadas 3D com Open Babel, e exporta PDB.

    Este metodo produz geometria razoavel para RNAs curtos mas
    sem considerar a estrutura secundaria.

    Args:
        sequence: Sequencia RNA/DNA.

    Returns:
        Conteudo PDB como string, ou None se Open Babel nao estiver
        disponivel.
    """
    try:
        from openbabel import openbabel as ob

        logger.info(
            "Construindo estrutura 3D com Open Babel para %d nt.",
            len(sequence),
        )

        # Gerar PDB usando coordenadas helicoidais e Open Babel
        # para otimizar a geometria
        pdb_lines = _build_helical_model(sequence)
        if not pdb_lines:
            return None

        # Usar Open Babel para otimizar geometria (force field)
        conv = ob.OBConversion()
        conv.SetInAndOutFormats("pdb", "pdb")

        mol = ob.OBMol()
        conv.ReadString(mol, pdb_lines)

        if mol.NumAtoms() == 0:
            logger.warning("Open Babel nao conseguiu ler o modelo PDB.")
            return None

        # Adicionar hidrogenios e otimizar
        mol.AddHydrogens()

        # Otimizacao com MMFF94 (se disponivel) ou UFF
        ff = ob.OBForceField.FindForceField("UFF")
        if ff is not None:
            ff.Setup(mol)
            ff.SteepestDescent(100)
            ff.GetCoordinates(mol)
            logger.info("Geometria otimizada com UFF force field.")

        # Exportar PDB otimizado
        pdb_optimized = conv.WriteString(mol)
        logger.info(
            "Open Babel gerou PDB com %d atomos.", mol.NumAtoms(),
        )
        return pdb_optimized

    except ImportError:
        logger.warning(
            "Open Babel nao disponivel. Tentando modelo helicoidal puro."
        )
        return None
    except Exception as exc:
        logger.warning(
            "Erro no Open Babel: %s. Usando modelo helicoidal.", exc,
        )
        return None


# ---------------------------------------------------------------------------
# Metodo 3: Modelo helicoidal A-form (fallback puro Python)
# ---------------------------------------------------------------------------


def _build_helical_model(sequence: str) -> str:
    """Constroi modelo 3D simplificado de RNA em conformacao A-form.

    Posiciona cada nucleotideo ao longo de uma helice com parametros
    da forma A do RNA. Cada nucleotideo e representado por 3 atomos
    pesados: P (fosfato), C4' (acucar), e N1/N9 (base).

    Este modelo NAO e preciso para docking de alta qualidade, mas
    permite:
    - Definir a grid box para docking
    - Identificar regioes de major/minor groove
    - Testar o pipeline completo sem dependencias externas

    Args:
        sequence: Sequencia RNA/DNA.

    Returns:
        String com conteudo PDB.
    """
    rna_seq = sequence.upper().replace("T", "U")
    lines: list[str] = []
    lines.append("REMARK   Generated by Marley - SL RNA helical model")
    lines.append(f"REMARK   Sequence: {rna_seq}")
    lines.append(f"REMARK   Length: {len(rna_seq)} nt")
    lines.append(
        f"REMARK   Date: {datetime.now(tz=timezone.utc).isoformat()}"
    )

    atom_serial = 1
    residue_num = 1

    for i, base in enumerate(rna_seq):
        # Angulo de rotacao para esta base
        theta = math.radians(TWIST_PER_BASE * i)
        # Altura ao longo do eixo da helice
        z = RISE_PER_BASE * i

        # Posicao do fosfato (P) - no backbone
        p_x = HELIX_RADIUS * math.cos(theta)
        p_y = HELIX_RADIUS * math.sin(theta)
        p_z = z

        # Posicao do acucar (C4') - ligeiramente para dentro
        sugar_radius = HELIX_RADIUS - 2.0
        c4_x = sugar_radius * math.cos(theta + math.radians(10))
        c4_y = sugar_radius * math.sin(theta + math.radians(10))
        c4_z = z + 0.5

        # Posicao da base (N1/N9) - projetada para fora
        base_radius = HELIX_RADIUS + 3.0
        base_atom = RNA_ATOM_TYPES.get(base, "N1")
        n_x = base_radius * math.cos(theta - math.radians(15))
        n_y = base_radius * math.sin(theta - math.radians(15))
        n_z = z + 1.0

        # Nome do residuo (3 letras)
        res_name = f"  {base}" if len(base) == 1 else base[:3]

        # Escrever ATOM records no formato PDB
        # Fosfato
        lines.append(
            f"ATOM  {atom_serial:5d}  P   {res_name:>3s} A{residue_num:4d}    "
            f"{p_x:8.3f}{p_y:8.3f}{p_z:8.3f}  1.00  0.00           P"
        )
        atom_serial += 1

        # Acucar C4'
        lines.append(
            f"ATOM  {atom_serial:5d}  C4' {res_name:>3s} A{residue_num:4d}    "
            f"{c4_x:8.3f}{c4_y:8.3f}{c4_z:8.3f}  1.00  0.00           C"
        )
        atom_serial += 1

        # Base nitrogenada
        lines.append(
            f"ATOM  {atom_serial:5d}  {base_atom:<3s} {res_name:>3s} A{residue_num:4d}    "
            f"{n_x:8.3f}{n_y:8.3f}{n_z:8.3f}  1.00  0.00           N"
        )
        atom_serial += 1

        residue_num += 1

    lines.append("END")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Conversao PDB -> PDBQT para receptor RNA
# ---------------------------------------------------------------------------


def _convert_pdb_to_pdbqt(pdb_path: Path, pdbqt_path: Path) -> bool:
    """Converte PDB do RNA para PDBQT (formato receptor AutoDock Vina).

    Para RNA como receptor, o PDBQT precisa de cargas parciais e
    tipos de atomo AutoDock. Usa Open Babel para a conversao.

    NOTA: RNA como receptor e menos comum que proteinas. O AutoDock Vina
    aceita PDBQT de RNA, mas a parametrizacao de cargas e aproximada.
    Para publicacao, recomenda-se usar rDock ou DOCK6 que tem
    parametros especificos para RNA.

    Args:
        pdb_path: Caminho do arquivo PDB de entrada.
        pdbqt_path: Caminho do arquivo PDBQT de saida.

    Returns:
        True se a conversao foi bem-sucedida.
    """
    try:
        from openbabel import openbabel as ob

        conv = ob.OBConversion()
        conv.SetInAndOutFormats("pdb", "pdbqt")
        # Opcoes para receptor (rigido, sem torsoes)
        conv.AddOption("r", ob.OBConversion.OUTOPTIONS)
        # Adicionar cargas Gasteiger
        conv.AddOption("p", ob.OBConversion.OUTOPTIONS)

        mol = ob.OBMol()
        conv.ReadFile(mol, str(pdb_path))

        if mol.NumAtoms() == 0:
            logger.error("PDB vazio: %s", pdb_path)
            return False

        # Adicionar hidrogenios polares (necessario para Vina)
        mol.AddHydrogens(True)  # polar only

        # Calcular cargas Gasteiger
        charge_model = ob.OBChargeModel.FindType("gasteiger")
        if charge_model is not None:
            charge_model.ComputeCharges(mol)

        conv.WriteFile(mol, str(pdbqt_path))
        logger.info(
            "Convertido PDB -> PDBQT: %s (%d atomos)",
            pdbqt_path,
            mol.NumAtoms(),
        )
        return True

    except ImportError:
        logger.warning(
            "Open Babel nao disponivel para conversao PDB->PDBQT. "
            "Gerando PDBQT simplificado."
        )
        return _fallback_pdb_to_pdbqt(pdb_path, pdbqt_path)


def _fallback_pdb_to_pdbqt(pdb_path: Path, pdbqt_path: Path) -> bool:
    """Conversao simplificada PDB -> PDBQT sem Open Babel.

    Adiciona colunas de carga parcial (0.000) e tipo AutoDock
    a cada ATOM record. Suficiente para testes, nao para producao.

    Args:
        pdb_path: Caminho do PDB de entrada.
        pdbqt_path: Caminho do PDBQT de saida.

    Returns:
        True se a conversao foi bem-sucedida.
    """
    # Mapeamento simplificado de elementos para tipos AutoDock
    ad_types: dict[str, str] = {
        "P": "P",
        "C": "C",
        "N": "NA",  # Nitrogenio aceptor de H
        "O": "OA",  # Oxigenio aceptor de H
        "S": "SA",
        "H": "HD",
    }

    try:
        with open(pdb_path, encoding="utf-8") as fh:
            pdb_lines = fh.readlines()

        pdbqt_lines: list[str] = []
        for line in pdb_lines:
            if line.startswith(("ATOM", "HETATM")):
                # Formato PDB: colunas 77-78 = elemento
                element = line[76:78].strip() if len(line) > 76 else "C"
                if not element:
                    # Tentar extrair do nome do atomo
                    atom_name = line[12:16].strip()
                    element = atom_name[0] if atom_name else "C"

                ad_type = ad_types.get(element, "C")
                charge = 0.000

                # PDBQT: mesma linha PDB + carga parcial + tipo AD
                pdbqt_line = (
                    f"{line[:54]}"
                    f"{charge:6.3f}"
                    f"    "
                    f"{ad_type:>2s}"
                )
                pdbqt_lines.append(pdbqt_line)
            elif line.startswith(("REMARK", "END")):
                pdbqt_lines.append(line.rstrip())

        pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
        with open(pdbqt_path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(pdbqt_lines))

        logger.info(
            "PDBQT simplificado gerado: %s (%d linhas ATOM)",
            pdbqt_path,
            sum(1 for l in pdbqt_lines if l.startswith("ATOM")),
        )
        return True

    except Exception as exc:
        logger.error("Falha na conversao PDB->PDBQT: %s", exc)
        return False


# ---------------------------------------------------------------------------
# Analise de binding pockets no RNA
# ---------------------------------------------------------------------------


def _identify_binding_pockets(
    pdb_content: str, dot_bracket: str, sequence: str,
) -> list[dict[str, Any]]:
    """Identifica potenciais bolsos de ligacao no RNA 3D.

    Em RNA, os bolsos de ligacao ocorrem em:
    1. Major groove (sulco maior) - acessivel a moleculas pequenas
    2. Minor groove (sulco menor) - mais estreito
    3. Bulges (protuberancias) - regioes nao pareadas
    4. Internal loops - regioes flexiveis
    5. Juncoes - onde stems se encontram

    Para o SL RNA (39nt, GC=28.2%), a estrutura e majoritariamente
    aberta, entao os bolsos estao nas regioes de loop/bulge.

    Args:
        pdb_content: Conteudo do arquivo PDB.
        dot_bracket: Estrutura secundaria.
        sequence: Sequencia do RNA.

    Returns:
        Lista de dicionarios com informacoes dos bolsos.
    """
    pockets: list[dict[str, Any]] = []

    # Extrair coordenadas dos atomos
    atoms: list[dict[str, float]] = []
    for line in pdb_content.split("\n"):
        if line.startswith("ATOM"):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append({"x": x, "y": y, "z": z})
            except (ValueError, IndexError):
                continue

    if not atoms:
        logger.warning("Nenhum atomo encontrado no PDB para analise de pockets.")
        return pockets

    # Calcular centroide geral
    cx = sum(a["x"] for a in atoms) / len(atoms)
    cy = sum(a["y"] for a in atoms) / len(atoms)
    cz = sum(a["z"] for a in atoms) / len(atoms)

    # Pocket 1: Regiao 5' do SL RNA (primeiros 13 nt)
    # Esta e a regiao mais conservada e funcional
    n_5prime = min(len(atoms), 39)  # 3 atomos por residuo, ~13 residuos
    if n_5prime > 0:
        p1_x = sum(a["x"] for a in atoms[:n_5prime]) / n_5prime
        p1_y = sum(a["y"] for a in atoms[:n_5prime]) / n_5prime
        p1_z = sum(a["z"] for a in atoms[:n_5prime]) / n_5prime
        pockets.append({
            "name": "5prime_exon",
            "description": "Regiao 5' do SL RNA (exon conservado, nt 1-13)",
            "center_x": round(p1_x, 1),
            "center_y": round(p1_y, 1),
            "center_z": round(p1_z, 1),
            "size_x": 25.0,
            "size_y": 25.0,
            "size_z": 25.0,
            "rationale": (
                "Regiao do mini-exon (nt 1-13) adicionada a todos os mRNAs "
                "por trans-splicing. Bloqueio aqui impede maturacao do mRNA."
            ),
        })

    # Pocket 2: Regiao central (splice site, nt 14-26)
    n_mid_start = min(39, len(atoms))  # atomo ~14*3
    n_mid_end = min(78, len(atoms))    # atomo ~26*3
    mid_atoms = atoms[n_mid_start:n_mid_end] if n_mid_end > n_mid_start else atoms
    if mid_atoms:
        p2_x = sum(a["x"] for a in mid_atoms) / len(mid_atoms)
        p2_y = sum(a["y"] for a in mid_atoms) / len(mid_atoms)
        p2_z = sum(a["z"] for a in mid_atoms) / len(mid_atoms)
        pockets.append({
            "name": "splice_site",
            "description": "Regiao do splice site (nt 14-26)",
            "center_x": round(p2_x, 1),
            "center_y": round(p2_y, 1),
            "center_z": round(p2_z, 1),
            "size_x": 25.0,
            "size_y": 25.0,
            "size_z": 25.0,
            "rationale": (
                "Sitio de splicing onde o SL RNA e clivado e transferido. "
                "Bloqueio aqui impede a reacao de trans-splicing."
            ),
        })

    # Pocket 3: Blind docking (centroide geral com box grande)
    pockets.append({
        "name": "blind_full",
        "description": "Docking cego cobrindo todo o RNA",
        "center_x": round(cx, 1),
        "center_y": round(cy, 1),
        "center_z": round(cz, 1),
        "size_x": 50.0,
        "size_y": 50.0,
        "size_z": 50.0,
        "rationale": (
            "Varredura completa para descobrir sitios de ligacao "
            "nao previstos pela analise de estrutura secundaria."
        ),
    })

    # Identificar regioes nao pareadas (potenciais bulges/loops)
    unpaired_runs: list[tuple[int, int]] = []
    i = 0
    while i < len(dot_bracket):
        if dot_bracket[i] == ".":
            start = i
            while i < len(dot_bracket) and dot_bracket[i] == ".":
                i += 1
            if (i - start) >= 3:  # Runs de 3+ bases abertas
                unpaired_runs.append((start, i - 1))
        else:
            i += 1

    for run_idx, (run_start, run_end) in enumerate(unpaired_runs[:3]):
        # Mapear residuos para atomos (3 atomos por residuo)
        atom_start = run_start * 3
        atom_end = min((run_end + 1) * 3, len(atoms))
        run_atoms = atoms[atom_start:atom_end]
        if run_atoms:
            pr_x = sum(a["x"] for a in run_atoms) / len(run_atoms)
            pr_y = sum(a["y"] for a in run_atoms) / len(run_atoms)
            pr_z = sum(a["z"] for a in run_atoms) / len(run_atoms)
            pockets.append({
                "name": f"loop_{run_idx + 1}",
                "description": (
                    f"Loop/bulge nao pareado (nt {run_start + 1}-{run_end + 1})"
                ),
                "center_x": round(pr_x, 1),
                "center_y": round(pr_y, 1),
                "center_z": round(pr_z, 1),
                "size_x": 20.0,
                "size_y": 20.0,
                "size_z": 20.0,
                "rationale": (
                    f"Regiao de loop com {run_end - run_start + 1} nucleotideos "
                    "nao pareados. Regioes flexiveis sao mais acessiveis "
                    "a moleculas pequenas."
                ),
            })

    logger.info("Identificados %d bolsos de ligacao no SL RNA.", len(pockets))
    return pockets


# ---------------------------------------------------------------------------
# API publica
# ---------------------------------------------------------------------------


def generate_sl_rna_3d(
    force: bool = False,
    dry_run: bool = False,
    method: str = "auto",
) -> str:
    """Gera estrutura 3D do SL RNA e prepara para docking.

    Pipeline:
    1. Carrega estrutura secundaria (dot-bracket) dos resultados anteriores.
    2. Gera estrutura 3D usando o metodo selecionado.
    3. Converte PDB para PDBQT (formato receptor Vina).
    4. Identifica bolsos de ligacao (binding pockets).
    5. Salva metadados para uso pelos modulos de docking.

    Args:
        force: Regenerar mesmo se arquivos ja existem.
        dry_run: Modo teste - gera modelo helicoidal sem dependencias.
        method: Metodo de modelagem: "auto", "rnacomposer", "openbabel",
                ou "helix". Default "auto" tenta todos em ordem.

    Returns:
        Caminho do arquivo de metadados JSON.
    """
    if SL_METADATA_PATH.exists() and not force:
        logger.info(
            "Estrutura 3D do SL RNA ja existe em %s. Use --force para regenerar.",
            SL_METADATA_PATH,
        )
        return str(SL_METADATA_PATH)

    # Criar diretorios
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Passo 1: Carregar estrutura secundaria
    dot_bracket = _load_existing_structure()
    logger.info("Estrutura secundaria: %s", dot_bracket)

    # Passo 2: Gerar estrutura 3D
    pdb_content: str | None = None
    method_used = "none"

    if dry_run:
        method = "helix"
        logger.info("[DRY RUN] Usando modelo helicoidal simplificado.")

    if method in ("auto", "rnacomposer"):
        pdb_content = _try_rnacomposer(SL_SEQUENCE, dot_bracket)
        if pdb_content:
            method_used = "rnacomposer"

    if pdb_content is None and method in ("auto", "openbabel"):
        pdb_content = _try_openbabel_build(SL_SEQUENCE)
        if pdb_content:
            method_used = "openbabel"

    if pdb_content is None or method == "helix":
        pdb_content = _build_helical_model(SL_SEQUENCE)
        method_used = "helix"
        logger.info(
            "Modelo helicoidal A-form gerado para %d nucleotideos.",
            SL_LENGTH,
        )

    # Salvar PDB
    with open(SL_PDB_PATH, "w", encoding="utf-8") as fh:
        fh.write(pdb_content)
    logger.info("PDB salvo em %s (metodo: %s)", SL_PDB_PATH, method_used)

    # Passo 3: Converter para PDBQT
    pdbqt_ok = _convert_pdb_to_pdbqt(SL_PDB_PATH, SL_PDBQT_PATH)
    if not pdbqt_ok:
        logger.warning(
            "Conversao PDB->PDBQT falhou. O PDB ainda pode ser usado "
            "com ferramentas alternativas (rDock, DOCK6)."
        )

    # Passo 4: Identificar binding pockets
    pockets = _identify_binding_pockets(pdb_content, dot_bracket, SL_SEQUENCE)

    # Passo 5: Salvar metadados
    metadata: dict[str, Any] = {
        "generated_at": datetime.now(tz=timezone.utc).isoformat(),
        "target": "SL RNA (Spliced Leader RNA)",
        "organism": "Leishmania infantum",
        "sequence": SL_SEQUENCE,
        "sequence_rna": SL_SEQUENCE.replace("T", "U"),
        "length_nt": SL_LENGTH,
        "gc_content": SL_GC_CONTENT,
        "dot_bracket": dot_bracket,
        "modeling_method": method_used,
        "pdb_path": str(SL_PDB_PATH),
        "pdbqt_path": str(SL_PDBQT_PATH) if pdbqt_ok else None,
        "pdbqt_ready": pdbqt_ok,
        "binding_pockets": pockets,
        "quality_notes": {
            "rnacomposer": (
                "Alta qualidade - modelagem por fragmentos baseada no PDB. "
                "Requer internet."
            ),
            "openbabel": (
                "Media qualidade - geometria otimizada com UFF force field. "
                "Nao considera estrutura secundaria."
            ),
            "helix": (
                "Baixa qualidade - modelo helicoidal A-form simplificado. "
                "Suficiente para testes e definicao de grid box. "
                "NAO recomendado para publicacao."
            ),
        },
    }

    with open(SL_METADATA_PATH, "w", encoding="utf-8") as fh:
        json.dump(metadata, fh, indent=2, ensure_ascii=False)
    logger.info("Metadados salvos em %s", SL_METADATA_PATH)

    # Resumo
    logger.info("=" * 60)
    logger.info("RESUMO - Estrutura 3D do SL RNA")
    logger.info("  Sequencia: %s (%d nt)", SL_SEQUENCE, SL_LENGTH)
    logger.info("  Metodo: %s", method_used)
    logger.info("  PDB: %s", SL_PDB_PATH)
    logger.info("  PDBQT: %s", SL_PDBQT_PATH if pdbqt_ok else "NAO GERADO")
    logger.info("  Binding pockets: %d", len(pockets))
    for pocket in pockets:
        logger.info(
            "    - %s: centro=(%.1f, %.1f, %.1f)",
            pocket["name"],
            pocket["center_x"],
            pocket["center_y"],
            pocket["center_z"],
        )
    logger.info("=" * 60)

    return str(SL_METADATA_PATH)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Gerar estrutura 3D do SL RNA de L. infantum "
            "para docking de moleculas pequenas."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Regenerar mesmo se arquivos ja existem.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Modo teste: usa modelo helicoidal simplificado.",
    )
    parser.add_argument(
        "--method",
        choices=["auto", "rnacomposer", "openbabel", "helix"],
        default="auto",
        help="Metodo de modelagem 3D (default: auto).",
    )
    args = parser.parse_args()

    result = generate_sl_rna_3d(
        force=args.force,
        dry_run=args.dry_run,
        method=args.method,
    )
    logger.info("Completo: %s", result)
