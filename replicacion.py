"""Simulación de la replicación semiconservativa del ADN."""

COMPLEMENTO_ADN = {"A": "T", "T": "A", "G": "C", "C": "G"}
LONGITUD_CEBADOR = 5
LONGITUD_OKAZAKI = 10
MAX_LONGITUD = 60


def complementaria(secuencia: str):
    """
    Devuelve la cadena complementaria de ADN, base a base.

    La cadena se escribe en el mismo orden que la original, sin invertirla, si la entrada se lee en sentido 5'->3',
    el resultado se lee 3'->5'.

    Args:
        - secuencia (str): Secuencia de ADN (solo A, T, G, C).

    Returns:
        - str: Secuencia complementaria (A-T, G-C) de la misma longitud.
    """
    return "".join(COMPLEMENTO_ADN[b] for b in secuencia)


def cebador_arn(secuencia_adn: str):
    """
    Convierte una secuencia de ADN en un cebador de ARN.

    Cambia cada T por U y escribe el resultado en minúsculas, que es la convención del simulador para distinguir el
    ARN del ADN. Por ejemplo, "TTGAC" se convierte en "uugac".

    Args:
        - secuencia_adn (str): Fragmento de ADN en mayúsculas.

    Returns:
        - str: Cebador de ARN en minúsculas.
    """
    return secuencia_adn.replace("T", "U").lower()


def sintetizar_rezagada(codificante: str):
    """
    Genera los fragmentos de Okazaki de la cadena rezagada.

    La hebra codificante actúa como molde de la rezagada. Se divide en tramos de LONGITUD_CEBADOR + LONGITUD_OKAZAKI
    nucleótidos (el último puede ser más corto) y cada tramo se complementa. El cebador de ARN ocupa los últimos
    LONGITUD_CEBADOR nucleótidos del fragmento, que es su extremo 5' (si el fragmento es más corto, el cebador es el
    fragmento entero).

    Args:
        - codificante (str): Hebra codificante de ADN (5'->3').

    Returns:
        - list[tuple]: Una tupla (inicio, fin, cebador, fragmento) por cada fragmento. Las posiciones cuentan desde 1
        sobre la codificante, el cebador es ARN en minúsculas y el fragmento es ADN escrito en sentido 3'->5'.
    """
    paso = LONGITUD_CEBADOR + LONGITUD_OKAZAKI
    fragmentos = []

    for inicio in range(0, len(codificante), paso):
        fin = min(inicio + paso, len(codificante))
        fragmento = complementaria(codificante[inicio:fin])
        cebador = cebador_arn(fragmento[-LONGITUD_CEBADOR:])
        fragmentos.append((inicio + 1, fin, cebador, fragmento))

    return fragmentos


def replicar_adn(codificante: str):
    """
    Replica la molécula y devuelve las dos moléculas hijas.

    La replicación es semiconservativa, cada hija conserva una hebra parental y recibe una hebra nueva. La cadena líder
    se obtiene complementando el molde y la rezagada se une a partir de los fragmentos de Okazaki.

    Args:
        - codificante (str): Hebra codificante de ADN (5'->3').

    Returns:
        - ((str, str), (str, str)): Las dos moléculas hijas, cada una como (hebra de arriba 5'->3', hebra de abajo 3'->5').
          Hija 1 = (líder nueva, molde parental).
          Hija 2 = (codificante parental, rezagada nueva).
    """
    molde = complementaria(codificante)
    lider = complementaria(molde)
    rezagada = "".join(
        fragmento
        for _, _, _, fragmento in sintetizar_rezagada(codificante)
    )

    return (lider, molde), (codificante, rezagada)


def mostrar_replicacion(codificante: str):
    """
    Muestra de forma resumida la replicación del ADN.

    Solo se dibujan los primeros MAX_LONGITUD nucleótidos, pero los cálculos usan la secuencia completa.

    Args:
        - codificante (str): Hebra codificante de ADN (5'->3').

    Returns:
        - ((str, str), (str, str)): Las dos moléculas hijas, igual que
          replicar_adn.
    """
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
