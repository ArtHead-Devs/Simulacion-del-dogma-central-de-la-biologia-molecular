"""
replicacion.py
--------------
Simula el Paso 3: Replicación del ADN (modelo semiconservativo).
Recibe una cadena codificante (5'->3') y muestra el proceso completo en consola.
"""

COMPLEMENTO_ADN = {"A": "T", "T": "A", "G": "C", "C": "G"}
LONGITUD_CEBADOR = 5
LONGITUD_OKAZAKI = 10


def complementaria(sec):
    """Devuelve la cadena complementaria de ADN (misma dirección)."""
    return "".join(COMPLEMENTO_ADN[b] for b in sec)


def replicar_adn(cadena_codificante):
    """
    Simula la replicación del ADN.
    Recibe la cadena codificante (5'->3') y devuelve (hija1, hija2).
    """
    cod = cadena_codificante.upper()
    n = len(cod)
    molde = complementaria(cod)

    print("\n=== REPLICACION DEL ADN ===\n")

    print("Molecula parental:")
    print(f"  Codificante  5' {cod} 3'")
    print(f"               {'|' * n}")
    print(f"  Molde        3' {molde} 5'")

    print("\nEnzimas:")
    print("  Topoisomerasa : relaja la tension antes de la apertura")
    print("  Helicasa      : separa las dos hebras")
    print("  SSB           : estabiliza las hebras simples")
    print("  Primasa       : sintetiza cebadores de ARN")
    print("  ADN pol III   : extiende los cebadores sintetizando ADN")
    print("  ADN pol I     : elimina cebadores y los reemplaza por ADN")
    print("  Ligasa        : une los fragmentos de Okazaki")

    mitad = n // 2
    print(f"\nApertura de la horquilla en posicion {mitad}:")
    print(f"  5' {cod[:mitad]}{' ' * (n - mitad)} 3'")
    print(f"     {'|' * mitad}")
    print(f"  3' {' ' * mitad}{molde[mitad:]} 5'")

    cebador_l = molde[:LONGITUD_CEBADOR].lower()
    ext_l = cod[LONGITUD_CEBADOR:]
    print("\nCadena LIDER (sintesis continua):")
    print(f"  Molde        3' {molde} 5'")
    print(f"  Cebador ARN  5' {cebador_l} 3'   (Primasa, {LONGITUD_CEBADOR} nt)")
    print(f"  Extension        {' ' * LONGITUD_CEBADOR}{ext_l}    (ADN pol III)")
    print(f"  Sin cebador  5' {cod} 3'   (ADN pol I elimina el cebador)")

    print("\nCadena REZAGADA (sintesis discontinua, fragmentos de Okazaki):")
    bloque = LONGITUD_CEBADOR + LONGITUD_OKAZAKI
    fragmentos = []
    pos, num = 0, 1
    while pos < n:
        fin = min(pos + bloque, n)
        seg = cod[pos:fin]
        nueva = complementaria(seg)
        cab_r = min(LONGITUD_CEBADOR, len(seg))
        cebador_r = nueva[:cab_r].lower()
        extension_r = nueva[cab_r:].upper()
        print(f"  Fragmento {num}  [{pos+1}..{fin}]")
        print(f"    Molde      3' {seg} 5'")
        print(f"    Cebador    5' {cebador_r} 3'   (Primasa)")
        if extension_r:
            print(f"    Extension     {' ' * cab_r}{extension_r}    (ADN pol III)")
        print(f"    Sin cebador5' {nueva.upper()} 3'   (ADN pol I)")
        fragmentos.append(nueva.upper())
        pos += bloque
        num += 1

    rezagada = "".join(fragmentos)
    print(f"\n  Ligasa une los fragmentos:")
    print(f"    5' {rezagada} 3'")

    hija1 = cod
    hija2 = cod

    print("\nMoleculas hijas (replicacion semiconservativa):")
    print(f"  Hija 1   5' {hija1} 3'  (hebra parental)")
    print(f"               {'|' * n}")
    print(f"           3' {molde} 5'  (hebra nueva)")
    print(f"  Hija 2   5' {hija2} 3'  (hebra nueva)")
    print(f"               {'|' * n}")
    print(f"           3' {complementaria(hija2)} 5'  (hebra parental)")

    print("\nReplicacion semiconservativa completada.\n")

    return hija1, hija2
