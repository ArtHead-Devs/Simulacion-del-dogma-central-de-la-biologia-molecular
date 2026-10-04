"""Simulación de la replicación semiconservativa del ADN."""

COMPLEMENTO_ADN = {"A": "T", "T": "A", "G": "C", "C": "G"}
LONGITUD_CEBADOR = 5
LONGITUD_OKAZAKI = 10
VISTA = 60
SANGRIA = 22


def complementaria(secuencia):
    """Devuelve la cadena complementaria de ADN (misma dirección)."""
    return "".join(COMPLEMENTO_ADN[b] for b in secuencia)


def cebador_arn(secuencia_adn):
    """Convierte ADN en un cebador de ARN (T pasa a U, en minúsculas)."""
    return secuencia_adn.replace("T", "U").lower()


def sintetizar_rezagada(codificante):
    """Sintetiza la cadena rezagada en fragmentos de Okazaki.

    Devuelve una lista de (inicio, fin, cebador, fragmento), con las
    posiciones contadas desde 1.
    """
    paso = LONGITUD_CEBADOR + LONGITUD_OKAZAKI
    fragmentos = []
    for inicio in range(0, len(codificante), paso):
        fin = min(inicio + paso, len(codificante))
        fragmento = complementaria(codificante[inicio:fin])
        cebador = cebador_arn(fragmento[-LONGITUD_CEBADOR:])
        fragmentos.append((inicio + 1, fin, cebador, fragmento))
    return fragmentos


def replicar_adn(codificante):
    """Replica la molécula y devuelve sus dos moléculas hijas.

    Cada hija es una tupla (hebra codificante, hebra molde).
    """
    molde = complementaria(codificante)
    lider = complementaria(molde)
    rezagada = "".join(f[3] for f in sintetizar_rezagada(codificante))
    return (lider, molde), (codificante, rezagada)


def fila(nombre, extremo_5, texto, extremo_3):
    """Imprime una hebra con sus extremos, recortada al ancho de vista."""
    puntos = " ..." if len(texto) > VISTA else ""
    print(f"  {nombre:<16} {extremo_5} {texto[:VISTA]}{puntos} {extremo_3}")


def barras(longitud, desde=0):
    """Imprime los puentes de hidrógeno entre dos hebras."""
    ancho = min(longitud, VISTA) - desde
    print(" " * (SANGRIA + desde) + "|" * ancho)


def mostrar_parental(codificante, molde):
    """Muestra la molécula de ADN de partida."""
    print("\n=== REPLICACIÓN DEL ADN ===\n")
    print("Molécula parental:")
    fila("Codificante", "5'", codificante, "3'")
    barras(len(codificante))
    fila("Molde", "3'", molde, "5'")


def mostrar_apertura(codificante, molde):
    """Muestra la apertura de la horquilla de replicación."""
    posicion = min(len(codificante), VISTA) // 2
    print(f"\nApertura [Topoisomerasa, Helicasa, SSB] en {posicion}:")
    fila("Codificante", "5'", codificante, "3'")
    barras(len(codificante), posicion)
    fila("Molde", "3'", molde, "5'")
    print(" " * (SANGRIA + posicion) + "^")


def mostrar_lider(molde, lider):
    """Muestra la síntesis continua de la cadena líder."""
    cebador = cebador_arn(lider[:LONGITUD_CEBADOR])
    print("\nCadena LÍDER [Primasa, ADN polimerasa III]:")
    print(f"  Cebador de ARN: {cebador}")
    fila("Molde", "3'", molde, "5'")
    barras(len(molde))
    fila("Líder (nueva)", "5'", cebador + lider[LONGITUD_CEBADOR:], "3'")


def rezagada_con_cebadores(fragmentos):
    """Une los fragmentos dejando los cebadores de ARN en minúsculas."""
    return "".join(
        fragmento[:len(fragmento) - len(cebador)] + cebador
        for _, _, cebador, fragmento in fragmentos
    )


def mostrar_rezagada(codificante, fragmentos):
    """Muestra la síntesis discontinua de la cadena rezagada."""
    limites = [" "] * len(codificante)
    for inicio, fin, _, _ in fragmentos:
        limites[fin - 1] = "]"
        limites[inicio - 1] = "["

    print("\nCadena REZAGADA [Primasa, ADN polimerasa III]:")
    fila("Molde", "5'", codificante, "3'")
    barras(len(codificante))
    fila("Rezagada (nueva)", "3'", rezagada_con_cebadores(fragmentos), "5'")
    print(f"  {'Fragmentos':<16}    {''.join(limites)[:VISTA]}\n")
    for numero, (inicio, fin, cebador, _) in enumerate(fragmentos[:3], 1):
        print(f"  Fragmento {numero} [{inicio}..{fin}] cebador: {cebador}")
    print(f"  Fragmentos de Okazaki: {len(fragmentos)}")
    print("  Cebadores a ADN [ADN polimerasa I], muescas [Ligasa]")


def mostrar_hijas(hija_1, hija_2):
    """Muestra las dos moléculas hijas."""
    lider, molde = hija_1
    codificante, rezagada = hija_2
    print("\nHija 1 = hebra parental (molde) + hebra LÍDER nueva")
    fila("Nueva (líder)", "5'", lider, "3'")
    barras(len(lider))
    fila("Parental", "3'", molde, "5'")
    print("\nHija 2 = hebra parental (codificante) + hebra REZAGADA nueva")
    fila("Parental", "5'", codificante, "3'")
    barras(len(codificante))
    fila("Nueva (rezagada)", "3'", rezagada, "5'")


def mostrar_replicacion(codificante):
    """Muestra todo el proceso de replicación y devuelve las hijas."""
    molde = complementaria(codificante)
    hija_1, hija_2 = replicar_adn(codificante)
    mostrar_parental(codificante, molde)
    mostrar_apertura(codificante, molde)
    mostrar_lider(molde, hija_1[0])
    mostrar_rezagada(codificante, sintetizar_rezagada(codificante))
    mostrar_hijas(hija_1, hija_2)
    return hija_1, hija_2