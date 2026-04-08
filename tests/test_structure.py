"""Tests for module 07 -- 3D structure prediction and visualization.

Tests cover the ESMFold API client, region mapping, pLDDT extraction,
and PyMOL / ChimeraX visualization script generation implemented in
``core.structure``.
"""

from __future__ import annotations

import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from core.structure import (
    COLOR_SCHEME,
    ConstructRegion,
    extract_plddt_scores,
    generate_chimerax_script,
    generate_pymol_script,
    map_construct_regions,
    predict_structure,
)

# ---------------------------------------------------------------------------
# Fixtures / constants
# ---------------------------------------------------------------------------

# Mock PDB with 5 ATOM records; 2 are CA atoms (residues 1 and 2).
MOCK_PDB = textwrap.dedent("""\
    ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00 89.50
    ATOM      2  CA  MET A   1      26.266  25.413   2.842  1.00 89.50
    ATOM      3  C   MET A   1      26.913  26.639   3.531  1.00 89.50
    ATOM      4  N   ASP A   2      27.886  26.463   4.417  1.00 91.20
    ATOM      5  CA  ASP A   2      28.621  27.524   5.068  1.00 91.20
    END
""")

MOCK_PDB_NO_ATOMS = textwrap.dedent("""\
    HEADER    MOCK STRUCTURE
    REMARK   THIS IS NOT A REAL PDB
    END
""")

# Construct component sequences (from pipeline.06_construct).
TPA = "MDAMKRGLCCVLLLCGAVFVSAS"  # 23 aa
L7L12 = (
    "MAKLSTDELLDAFKEMTLLELSDFVKKFEETFEVTAAAPVAVAAAGAAPAGAAVEAAEEQ"
    "SEFDVILEAAGDKKIGVIKVVREIVSGLGLKEAKDLVDGAPKPLLEKVAKEAADEAKAKL"
    "EAAGATVTVK"
)  # 130 aa
EAAAK = "EAAAK"  # 5 aa
AAY = "AAY"
GPGPG = "GPGPG"


def _build_construct(
    ctl_epitopes: list[str],
    htl_epitopes: list[str] | None = None,
) -> str:
    """Build a construct string following the known assembly order."""
    parts = [TPA, L7L12, EAAAK]
    parts.append(AAY.join(ctl_epitopes))
    if htl_epitopes:
        parts.append(GPGPG)
        parts.append(GPGPG.join(htl_epitopes))
    return "".join(parts)


def _epitope_dict(sequence: str, epitope_type: str = "CTL") -> dict:
    """Build a minimal epitope dict matching the interface expected by
    ``map_construct_regions``."""
    return {"sequence": sequence, "epitope_type": epitope_type}


@pytest.fixture()
def sample_pdb_file(tmp_path: Path) -> Path:
    """Write a small mock PDB to a temporary file and return its path."""
    pdb_path = tmp_path / "structure.pdb"
    pdb_path.write_text(MOCK_PDB)
    return pdb_path


@pytest.fixture()
def sample_regions() -> list[ConstructRegion]:
    """A minimal list of ConstructRegion objects for script-generation tests."""
    return [
        ConstructRegion(
            name="signal_peptide",
            region_type="signal_peptide",
            start=1,
            end=23,
            color="#3B82F6",
            label="Signal peptide (tPA)",
        ),
        ConstructRegion(
            name="adjuvant",
            region_type="adjuvant",
            start=24,
            end=153,
            color="#22C55E",
            label="Adjuvant (L7L12)",
        ),
        ConstructRegion(
            name="CTL_1",
            region_type="ctl_epitope",
            start=159,
            end=169,
            color="#EF4444",
            label="CTL epitope 1",
        ),
    ]


# ---------------------------------------------------------------------------
# ESMFold client tests
# ---------------------------------------------------------------------------


class TestPredictStructure:
    """Tests for the ESMFold HTTP client wrapper."""

    @patch("core.structure.time.sleep")
    @patch("core.structure.requests.post")
    def test_predict_structure_success(
        self, mock_post: MagicMock, _mock_sleep: MagicMock, tmp_path: Path
    ) -> None:
        """A valid PDB response is written to the output file."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = MOCK_PDB
        mock_response.raise_for_status = MagicMock()
        mock_post.return_value = mock_response

        out_path = tmp_path / "predicted.pdb"
        result = predict_structure("MKFV", output_path=out_path)

        assert result is not None
        assert out_path.exists()
        contents = out_path.read_text()
        assert "ATOM" in contents
        mock_post.assert_called_once()

    @patch("core.structure.time.sleep")
    @patch("core.structure.requests.post")
    def test_predict_structure_timeout(
        self, mock_post: MagicMock, _mock_sleep: MagicMock, tmp_path: Path
    ) -> None:
        """A Timeout exception results in None, no crash."""
        import requests as _req

        mock_post.side_effect = _req.exceptions.Timeout("timed out")

        result = predict_structure("MKFV", output_path=tmp_path / "out.pdb")

        assert result is None

    @patch("core.structure.time.sleep")
    @patch("core.structure.requests.post")
    def test_predict_structure_server_error(
        self, mock_post: MagicMock, _mock_sleep: MagicMock, tmp_path: Path
    ) -> None:
        """A 500 server error results in None."""
        import requests as _req

        mock_response = MagicMock()
        mock_response.status_code = 500
        mock_response.text = "Internal Server Error"
        mock_response.raise_for_status.side_effect = _req.exceptions.HTTPError(
            "500 Server Error"
        )
        mock_post.return_value = mock_response

        result = predict_structure("MKFV", output_path=tmp_path / "out.pdb")

        assert result is None

    def test_predict_structure_sequence_too_long(self, tmp_path: Path) -> None:
        """Sequences longer than 800 aa are rejected without an API call."""
        long_seq = "M" * 801

        with patch("core.structure.requests.post") as mock_post:
            result = predict_structure(
                long_seq, output_path=tmp_path / "out.pdb"
            )

        assert result is None
        mock_post.assert_not_called()

    @patch("core.structure.time.sleep")
    @patch("core.structure.requests.post")
    def test_predict_structure_validates_pdb(
        self, mock_post: MagicMock, _mock_sleep: MagicMock, tmp_path: Path
    ) -> None:
        """A response with no ATOM records is treated as invalid."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = MOCK_PDB_NO_ATOMS
        mock_response.raise_for_status = MagicMock()
        mock_post.return_value = mock_response

        result = predict_structure("MKFV", output_path=tmp_path / "out.pdb")

        assert result is None


# ---------------------------------------------------------------------------
# Region mapper tests
# ---------------------------------------------------------------------------


class TestMapConstructRegions:
    """Tests for mapping construct protein regions to residue boundaries."""

    def test_map_construct_regions_basic(self) -> None:
        """Verify region boundaries for a construct with 2 CTL epitopes."""
        epi1 = "KLFPGDEIFSV"  # 11 aa
        epi2 = "YMLDIFHEV"   # 9 aa
        construct = _build_construct([epi1, epi2])

        epitopes = [
            _epitope_dict(epi1, "CTL"),
            _epitope_dict(epi2, "CTL"),
        ]
        regions = map_construct_regions(
            construct,
            signal_peptide_name="tPA",
            adjuvant_name="L7L12",
            epitopes=epitopes,
        )

        region_names = [r.name for r in regions]

        assert "signal_peptide" in region_names
        assert "adjuvant" in region_names
        assert "linker_EAAAK" in region_names

        # Signal peptide: residues 1..23
        sp = next(r for r in regions if r.name == "signal_peptide")
        assert sp.start == 1
        assert sp.end == len(TPA)

        # Adjuvant: starts right after signal peptide
        adj = next(r for r in regions if r.name == "adjuvant")
        assert adj.start == len(TPA) + 1
        assert adj.end == len(TPA) + len(L7L12)

        # EAAAK linker: 5 residues right after adjuvant
        lnk = next(r for r in regions if r.name == "linker_EAAAK")
        assert lnk.end - lnk.start + 1 == 5

        # Both CTL epitopes must appear as regions.
        ctl_regions = [r for r in regions if r.name.startswith("CTL_")]
        assert len(ctl_regions) == 2

    def test_map_construct_regions_with_htl(self) -> None:
        """Include HTL epitopes and verify GPGPG linkers and HTL regions."""
        ctl_epi = "KLFPGDEIFSV"
        htl1 = "ACDEFGHIKLMNPQR"  # 15 aa
        htl2 = "CDEFGHIKLMNPQRS"  # 15 aa
        construct = _build_construct(
            [ctl_epi], htl_epitopes=[htl1, htl2]
        )

        epitopes = [
            _epitope_dict(ctl_epi, "CTL"),
            _epitope_dict(htl1, "HTL"),
            _epitope_dict(htl2, "HTL"),
        ]
        regions = map_construct_regions(
            construct,
            signal_peptide_name="tPA",
            adjuvant_name="L7L12",
            epitopes=epitopes,
        )

        region_names = [r.name for r in regions]

        # GPGPG bridge between CTL and HTL blocks.
        assert any("GPGPG" in n for n in region_names)

        # HTL epitopes should be represented.
        htl_regions = [r for r in regions if r.name.startswith("HTL_")]
        assert len(htl_regions) == 2

    def test_map_construct_regions_colors(self) -> None:
        """Each region type must receive a color from COLOR_SCHEME."""
        epi1 = "KLFPGDEIFSV"
        construct = _build_construct([epi1])

        epitopes = [_epitope_dict(epi1, "CTL")]
        regions = map_construct_regions(
            construct,
            signal_peptide_name="tPA",
            adjuvant_name="L7L12",
            epitopes=epitopes,
        )

        for region in regions:
            assert region.color, (
                f"Region '{region.name}' has no color"
            )
            # The color should be one of the values in COLOR_SCHEME.
            assert region.color in COLOR_SCHEME.values(), (
                f"Unexpected color '{region.color}' for '{region.name}'"
            )


# ---------------------------------------------------------------------------
# pLDDT extraction tests
# ---------------------------------------------------------------------------


class TestExtractPlddtScores:
    """Tests for B-factor / pLDDT extraction from PDB ATOM records."""

    def test_extract_plddt_scores(self) -> None:
        """Known B-factors for CA atoms in the mock PDB are extracted."""
        scores = extract_plddt_scores(MOCK_PDB)

        # The mock has 2 CA atoms: residue 1 (89.50) and residue 2 (91.20).
        assert len(scores) == 2
        # Each entry is (residue_num, plddt).
        assert scores[0] == (1, pytest.approx(89.50))
        assert scores[1] == (2, pytest.approx(91.20))

    def test_extract_plddt_empty(self) -> None:
        """A PDB with no ATOM records yields an empty list."""
        scores = extract_plddt_scores(MOCK_PDB_NO_ATOMS)
        assert scores == []


# ---------------------------------------------------------------------------
# PyMOL script tests
# ---------------------------------------------------------------------------


class TestGeneratePymolScript:
    """Tests for PyMOL .pml visualization script generation."""

    def test_generate_pymol_script(
        self, sample_regions: list[ConstructRegion], tmp_path: Path
    ) -> None:
        """The .pml script contains select and color commands for each region."""
        out = tmp_path / "vis.pml"
        generate_pymol_script("model.pdb", sample_regions, out)
        script = out.read_text()

        for region in sample_regions:
            assert region.name in script
            assert "select" in script
            assert "color" in script

    def test_generate_pymol_script_loads_pdb(
        self, sample_regions: list[ConstructRegion], tmp_path: Path
    ) -> None:
        """The script must start with a load command for the PDB file."""
        out = tmp_path / "vis.pml"
        generate_pymol_script("model.pdb", sample_regions, out)
        script = out.read_text()

        lines = [ln.strip() for ln in script.splitlines() if ln.strip()]
        assert lines[0].startswith("load"), (
            "PyMOL script should begin with a 'load' command"
        )
        assert "model.pdb" in lines[0]


# ---------------------------------------------------------------------------
# ChimeraX script tests
# ---------------------------------------------------------------------------


class TestGenerateChimeraxScript:
    """Tests for ChimeraX .cxc visualization script generation."""

    def test_generate_chimerax_script(
        self, sample_regions: list[ConstructRegion], tmp_path: Path
    ) -> None:
        """The .cxc script contains color commands for each region."""
        out = tmp_path / "vis.cxc"
        generate_chimerax_script("model.pdb", sample_regions, out)
        script = out.read_text()

        # File should be loadable.
        assert "open" in script.lower() or "model.pdb" in script
        for region in sample_regions:
            # Residue range should appear.
            assert str(region.start) in script
            assert str(region.end) in script
        # At least one color command.
        assert "color" in script.lower()
