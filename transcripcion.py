"""
transcripcion.py
----------------
Simula la transcripcion de ADN a ARNm.
Recibe una cadena codificante (5'->3') y muestra el proceso en consola.
"""
COMPLEMENTO_ADN = {"A": "T", "T": "A", "G": "C", "C": "G"}
COMPLEMENTO_TRANSCRIPCION = {"A": "U", "T": "A", "C": "G", "G": "C"}
VISTA = 60


def complementaria_adn(sec):
    """Devuelve la cadena de ADN complementaria (misma direccion de escritura)."""
    return "".join(COMPLEMENTO_ADN[b] for b in sec)


def transcribir_adn(cadena_codificante):
    """
    Simula la transcripcion del ADN a ARNm.
    Recibe la cadena codificante (5'->3') y devuelve el ARNm (5'->3').
    """
    cod = cadena_codificante.upper()
    n = len(cod)
    molde = complementaria_adn(cod)
    arnm = "".join(COMPLEMENTO_TRANSCRIPCION[b] for b in molde)

    v = min(n, VISTA)
    puntos = " ..." if n > VISTA else ""

    print("\n=== TRANSCRIPCION DE ADN a ARNm ===\n")

    print("Enzimas y moleculas:")
    print("  Promotor: señal de inicio donde se une la ARN polimerasa")
    print("  ARN polimerasa: lee el molde 3'->5' y sintetiza el ARNm 5'->3'")
    print("  Ribonucleotidos: ATP, UTP, GTP, CTP")
    print("  Terminador: senal de fin de la transcripcion")

    print("\nComplementariedad molde -> ARNm:  A->U   T->A   C->G   G->C\n")

    if n > VISTA:
        print(f"(Se dibujan los primeros {VISTA} de {n} nt; el calculo usa la secuencia completa)\n")
    print(f"  Codificante  5' {cod[:v]}{puntos} 3'")
    print(f"                  {'|' * v}")
    print(f"  Molde        3' {molde[:v]}{puntos} 5'")
    print(f"                  {'|' * v}")
    print(f"  ARNm         5' {arnm[:v]}{puntos} 3'")

    return arnm