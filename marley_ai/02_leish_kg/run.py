"""Execucao do modulo 02_leish_kg — Knowledge Graph de Leishmania.

Uso:
    python -m marley_ai.02_leish_kg.run
"""

from __future__ import annotations

from typing import Any

from marley_ai.config import AIModuleConfig
from marley_ai.envelope import Timer, create_envelope, write_result


def main(config: AIModuleConfig | None = None) -> dict[str, Any]:
    """Constroi o Knowledge Graph de Leishmania.

    TODO: Implementar:
        1. Baixar proteoma de L. infantum JPCM5 do TriTrypDB
        2. Extrair GO terms e vias metabolicas do UniProt
        3. Integrar epitopos do vaccine_platforms/shared/epitopes.py
        4. Integrar alvos ASO do aso_math/results/
        5. Extrair relacoes de abstracts PubMed via PubMedBERT
        6. Construir grafo NetworkX com tipagem de nos e arestas
        7. Serializar em formato GraphML para visualizacao

    Args:
        config: Configuracao do modulo. Se None, usa defaults.

    Returns:
        Envelope com estatisticas do grafo.
    """
    envelope = create_envelope("02_leish_kg")

    with Timer() as timer:
        # TODO: Ingestao de dados do TriTrypDB
        # API REST: https://tritrypdb.org/tritrypdb/service/
        # Foco: L. infantum JPCM5 (referencia do projeto)

        # TODO: Anotacoes UniProt
        # Mapear gene_ids LINF_* para UniProt accessions
        # Extrair GO terms, domínios Pfam, vias KEGG

        # TODO: Integracao com pipeline Marley
        # Importar epitopos de vaccine_platforms.shared.epitopes
        # Importar resultados ASO de aso_math.results/

        # TODO: Extracao de relacoes via NLP
        # PubMed query: "Leishmania infantum" + termos de interesse
        # Modelo: PubMedBERT fine-tuned para relacoes biomedicas

        # TODO: Construcao do grafo
        # Tipos de nos: protein, gene, pathway, drug, epitope, disease
        # Tipos de arestas: encodes, participates_in, targets, binds, etc.

        envelope["status"] = "stub"
        envelope["summary"]["conclusion"] = (
            "Knowledge Graph em construcao. "
            "Integrara proteoma, vias metabolicas e literatura."
        )
        envelope["summary"]["key_metrics"] = {
            "n_nodes": 0,
            "n_edges": 0,
            "n_relation_types": 0,
        }

    envelope["runtime_seconds"] = timer.elapsed
    output_path = write_result(envelope)
    print(f"[02_leish_kg] Resultado salvo em {output_path}")

    return envelope


if __name__ == "__main__":
    main()
