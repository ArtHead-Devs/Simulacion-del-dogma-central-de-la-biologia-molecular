"""Simulación de la replicación semiconservativa del ADN."""

COMPLEMENTO_ADN = {"A": "T", "T": "A", "G": "C", "C": "G"}
LONGITUD_CEBADOR = 5
LONGITUD_OKAZAKI = 10
MAX_LONGITUD = 60


def complementaria(secuencia):
    """Devuelve la cadena complementaria de ADN."""
    return "".join(COMPLEMENTO_ADN[b] for b in secuencia)


def cebador_arn(secuencia_adn):
    """Convierte una secuencia de ADN en un cebador de ARN."""
    return secuencia_adn.replace("T", "U").lower()


def sintetizar_rezagada(codificante):
    """Genera los fragmentos de Okazaki de la cadena rezagada."""
    paso = LONGITUD_CEBADOR + LONGITUD_OKAZAKI
    fragmentos = []

    for inicio in range(0, len(codificante), paso):
        fin = min(inicio + paso, len(codificante))
        fragmento = complementaria(codificante[inicio:fin])
        cebador = cebador_arn(fragmento[-LONGITUD_CEBADOR:])
        fragmentos.append((inicio + 1, fin, cebador, fragmento))

    return fragmentos


def replicar_adn(codificante):
    """Replica la molécula y devuelve las dos moléculas hijas."""
    molde = complementaria(codificante)
    lider = complementaria(molde)
    rezagada = "".join(
        fragmento
        for _, _, _, fragmento in sintetizar_rezagada(codificante)
    )

    return (lider, molde), (codificante, rezagada)


def mostrar_replicacion(codificante):
    """Muestra de forma resumida la replicación del ADN."""
    molde = complementaria(codificante)
    hija_1, hija_2 = replicar_adn(codificante)
    fragmentos = sintetizar_rezagada(codificante)

    ancho = min(len(codificante), MAX_LONGITUD)
    puntos = " ..." if len(codificante) > MAX_LONGITUD else ""

    print("\n=== REPLICACIÓN DEL ADN ===")

    print("\n1. ADN parental:")
    print(f"  5' {codificante[:ancho]}{puntos} 3'")
    print(f"  3' {molde[:ancho]}{puntos} 5'")

    print("\n2. Apertura de la doble hélice [Topoisomerasa, Helicasa]:")
    print("  La topoisomerasa reduce la tensión y la helicasa separa las dos hebras.")
    print(f"  Hebra parental 1: 5' {codificante[:ancho]}{puntos} 3'")
    print(f"  Hebra parental 2: 3' {molde[:ancho]}{puntos} 5'")

    print("\n3. Cadena líder [Primasa, ADN polimerasa]:")
    cebador = cebador_arn(hija_1[0][:LONGITUD_CEBADOR])
    nueva_lider = cebador + hija_1[0][LONGITUD_CEBADOR:]
    print(f"  Cebador: {cebador}")
    print(f"  Nueva hebra: 5' {nueva_lider[:ancho]}{puntos} 3'")
    print("  Síntesis continua por ADN polimerasa.")

    print("\n4. Cadena rezagada [Primasa, ADN polimerasa]:")
    for numero, (inicio, fin, cebador, _) in enumerate(fragmentos[:3], 1):
        print(f"  Okazaki {numero} [{inicio}-{fin}]: cebador {cebador}")
    print(f"  Total de fragmentos: {len(fragmentos)}")
    print("  Síntesis discontinua por ADN polimerasa.")
    print("  La ADN polimerasa elimina los cebadores y los sustituye por ADN.")
    print("  La ligasa une los fragmentos de Okazaki.")

    print("\n5. Moléculas hijas:")

    print("  Hija 1:")
    print(f"    5' {hija_1[0][:ancho]}{puntos} 3'  [nueva, líder]")
    print(f"    3' {hija_1[1][:ancho]}{puntos} 5'  [parental]")

    print("  Hija 2:")
    print(f"    5' {hija_2[0][:ancho]}{puntos} 3'  [parental]")
    print(f"    3' {hija_2[1][:ancho]}{puntos} 5'  [nueva, rezagada]")

    return hija_1, hija_2
