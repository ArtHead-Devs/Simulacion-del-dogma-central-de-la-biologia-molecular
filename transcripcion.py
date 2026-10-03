"""
transcripcion.py
----------------
Simula el Paso 4: Transcripción de ADN a ARNm.
Recibe una cadena codificante (5'->3') y muestra el proceso completo en consola.
"""

COMPLEMENTO_ADN = {"A": "T", "T": "A", "G": "C", "C": "G"}
COMPLEMENTO_TRANSCRIPCION = {"A": "U", "T": "A", "C": "G", "G": "C"}


def complementaria_adn(sec):
    """Devuelve la cadena de ADN complementaria (misma dirección)."""
    return "".join(COMPLEMENTO_ADN[b] for b in sec)


def transcribir_adn(cadena_codificante):
    """
    Simula la transcripción del ADN a ARNm.
    Recibe la cadena codificante (5'->3') y devuelve el ARNm (5'->3').
    """
    cod = cadena_codificante.upper()
    molde = complementaria_adn(cod)
    arnm = "".join(COMPLEMENTO_TRANSCRIPCION[b] for b in molde)

    print("\n=== TRANSCRIPCION DE ADN a ARNm ===\n")

    print("Enzimas y moleculas:")
    print("  Promotor        : secuencia de inicio donde se une la ARN polimerasa")
    print("  ARN polimerasa  : lee la hebra molde 3' a 5' y sintetiza ARNm 5' a 3'")
    print("  Ribonucleotidos : ATP, UTP, GTP, CTP")
    print("  Terminador      : secuencia de fin de transcripcion")

    print("\nReglas de complementariedad (molde a ARNm):")
    print("  A del molde produce U en el ARNm")
    print("  T del molde produce A en el ARNm")
    print("  C del molde produce G en el ARNm")
    print("  G del molde produce C en el ARNm")

    print("\nSintesis del ARNm:")
    print(f"  Molde  3' {molde} 5'")
    print(f"  ARNm   5' {arnm} 3'")

    print("\nAlineamiento de las tres cadenas:")
    print(f"  Codificante  5' {cod}  3'")
    print(f"  Molde        3' {molde}  5'")
    print(f"  ARNm         5' {arnm}  3'")

    verificacion = cod.replace("T", "U")
    estado = "correcto" if arnm == verificacion else "error: revisar logica"
    print(f"\nVerificacion (ARNm igual a codificante sustituyendo T por U): {estado}\n")

    return arnm
