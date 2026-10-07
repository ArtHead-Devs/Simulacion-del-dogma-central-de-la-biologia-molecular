"""Simulación de la transcripción de ADN a ARN mensajero."""

from replicacion import MAX_LONGITUD, complementaria

COMPLEMENTO_TRANSCRIPCION = {"A": "U", "T": "A", "C": "G", "G": "C"}


def transcribir_adn(cadena_codificante: str):
    """
    Transcribe la hebra codificante de ADN a ARNm.

    Primero obtiene la hebra molde (complementaria de la codificante) y después genera el ARNm complementando el molde
    con las reglas A->U, T->A, C->G, G->C. El resultado tiene la misma secuencia que la codificante, con U en lugar de T.

    Args:
        - cadena_codificante (str): Hebra codificante de ADN (5'->3'). Se convierte a mayúsculas.

    Returns:
        - str: ARNm (5'->3').
    """
    codificante = cadena_codificante.upper()
    molde = complementaria(codificante)

    return "".join(COMPLEMENTO_TRANSCRIPCION[base] for base in molde)


def mostrar_transcripcion(codificante: str):
    """
    Muestra el proceso de transcripción.

    Imprime la regla de complementariedad y tres líneas, la hebra codificante, la hebra molde y el ARNm. Solo se
    imprimen los primeros MAX_LONGITUD nucleótidos, pero el cálculo usa la secuencia completa.

    Args:
        - codificante (str): Hebra codificante de ADN (5'->3'). Se convierte a mayúsculas.

    Returns:
        - str: ARNm (5'->3') completo.
    """
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
