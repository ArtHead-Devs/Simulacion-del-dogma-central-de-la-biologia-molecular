"""Simulación de la transcripción de ADN a ARN mensajero."""

from replicacion import MAX_LONGITUD, complementaria

COMPLEMENTO_TRANSCRIPCION = {"A": "U", "T": "A", "C": "G", "G": "C"}


def transcribir_adn(cadena_codificante):
    """Transcribe ADN a ARNm."""
    codificante = cadena_codificante.upper()
    molde = complementaria(codificante)

    return "".join(COMPLEMENTO_TRANSCRIPCION[base] for base in molde)


def mostrar_transcripcion(codificante):
    """Muestra el proceso de transcripción."""
    codificante = codificante.upper()
    molde = complementaria(codificante)
    arnm = transcribir_adn(codificante)

    ancho = min(len(codificante), MAX_LONGITUD)
    puntos = " ..." if len(codificante) > MAX_LONGITUD else ""

    print("\n=== TRANSCRIPCIÓN DE ADN A ARNm [ARN polimerasa] ===")
    print("  La ARN polimerasa lee la hebra molde (3'->5') y sintetiza el ARNm (5'->3').")
    print("  Complementariedad molde -> ARNm:  A->U  T->A  C->G  G->C\n")

    print(f"  ADN Codificante: 5' {codificante[:ancho]}{puntos} 3'")
    print(f"  ADN Molde:       3' {molde[:ancho]}{puntos} 5'")
    print(f"  ARNm:            5' {arnm[:ancho]}{puntos} 3'")

    return arnm