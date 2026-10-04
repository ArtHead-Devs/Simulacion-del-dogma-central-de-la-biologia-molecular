"""Simulación de la transcripción de ADN a ARN mensajero."""

from Bio.Seq import Seq

from replicacion import barras, complementaria, fila


def transcribir_adn(codificante):
    """Transcribe la hebra codificante (5'->3') a ARNm (5'->3').

    La ARN polimerasa lee el molde (3'->5') y sintetiza el ARNm (5'->3').
    """
    molde = complementaria(codificante)
    return str(Seq(molde).complement_rna())


def mostrar_transcripcion(codificante):
    """Muestra el proceso de transcripción y devuelve el ARNm."""
    molde = complementaria(codificante)
    arnm = transcribir_adn(codificante)
    print("\n=== TRANSCRIPCIÓN DE ADN A ARNm [ARN polimerasa] ===\n")
    print("Complementariedad molde -> ARNm:  A->U  T->A  C->G  G->C\n")
    fila("Codificante", "5'", codificante, "3'")
    barras(len(codificante))
    fila("Molde", "3'", molde, "5'")
    barras(len(molde))
    fila("ARNm", "5'", arnm, "3'")
    return arnm