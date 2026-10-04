"""
replicacion.py
--------------
Simula la replicacion del ADN (modelo semiconservativo).
Recibe una cadena codificante (5'->3') y muestra el proceso en consola.
La horquilla avanza hacia la derecha. MAYUSCULAS = ADN, minusculas = ARN (cebadores).
"""
from Bio.Seq import Seq

LONGITUD_CEBADOR = 5
LONGITUD_OKAZAKI = 10
VISTA = 60

def complementaria(sec):
    """Devuelve la cadena complementaria de ADN (misma direccion de escritura)."""
    return str(Seq(sec).complement())


def cebador_arn(sec_adn):
    """Un cebador es de ARN: la T pasa a U. Se escribe en minusculas."""
    return sec_adn.replace("T", "U").lower()


def replicar_adn(cadena_codificante):
    """Simula la replicacion. Recibe la cadena codificante (5'->3') y devuelve (hija1, hija2)."""
    cod = cadena_codificante.upper()
    n = len(cod)
    molde = complementaria(cod)
    v = min(n, VISTA)
    puntos = " ..." if n > VISTA else ""

    def fila(nombre, e5, texto, e3):
        print(f"  {nombre:<16} {e5} {texto[:v]}{puntos} {e3}")

    def barras(desde=0):
        print(" " * (22 + desde) + "|" * (v - desde))

    print("\n=== REPLICACION DEL ADN ===\n")
    if n > VISTA:
        print(f"(Se dibujan los primeros {VISTA} de {n} nt; el calculo usa la secuencia completa)\n")

    print("Molecula parental:")
    fila("Codificante", "5'", cod, "3'")
    barras()
    fila("Molde", "3'", molde, "5'")

    print("\nEnzimas:")
    print("  Topoisomerasa: relaja la tension antes de la apertura")
    print("  Helicasa: separa las dos hebras")
    print("  Primasa: sintetiza cebadores de ARN")
    print("  ADN polimerasa: extiende los cebadores (5'->3') y los sustituye por ADN")
    print("  Ligasa: une los fragmentos de Okazaki")

    h = v // 2
    print(f"\nApertura de la horquilla (posicion {h}, avanza hacia la derecha):")
    fila("Codificante", "5'", cod, "3'")
    barras(h)
    fila("Molde", "3'", molde, "5'")
    print(" " * (22 + h) + "^")

    cebador_l = cebador_arn(cod[:LONGITUD_CEBADOR])
    lider = cebador_l + cod[LONGITUD_CEBADOR:]
    print("\nCadena LIDER (sintesis continua, mismo sentido que la horquilla):")
    print(f"  Un unico cebador de ARN ({cebador_l}); la ADN polimerasa extiende sin parar.")
    fila("Molde", "3'", molde, "5'")
    barras()
    fila("Lider (nueva)", "5'", lider, "3'")

    print("\nCadena REZAGADA (sintesis discontinua, fragmentos de Okazaki):")
    print("  Molde = hebra codificante. Cada fragmento lleva su cebador en el extremo 5'")
    print("  (a la derecha) y se sintetiza hacia la izquierda, alejandose de la horquilla.")
    marcada, limites, fragmentos = "", [" "] * n, []
    for ini in range(0, n, LONGITUD_CEBADOR + LONGITUD_OKAZAKI):
        fin = min(ini + LONGITUD_CEBADOR + LONGITUD_OKAZAKI, n)
        nueva = complementaria(cod[ini:fin])
        cab = min(LONGITUD_CEBADOR, len(nueva))
        cebador = cebador_arn(nueva[-cab:])
        marcada += nueva[:-cab] + cebador
        limites[ini], limites[fin - 1] = "[", ("]" if fin - 1 > ini else "[")
        fragmentos.append((ini + 1, fin, cebador, nueva))
    fila("Molde", "5'", cod, "3'")
    barras()
    fila("Rezagada (nueva)", "3'", marcada, "5'")
    print(f"  {'Fragmentos':<16}    {''.join(limites)[:v]}\n")
    for num, (ini, fin, cebador, _) in enumerate(fragmentos[:3], 1):
        print(f"  Fragmento {num}  [{ini}..{fin}]  cebador de ARN: {cebador}")
    print(f"  Total: {len(fragmentos)} fragmentos de Okazaki (1 cebador por fragmento)")
    print("  La ADN polimerasa cambia los cebadores por ADN y la ligasa sella las muescas.")

    lider_final = lider.upper().replace("U", "T")
    rezagada_final = "".join(f[3] for f in fragmentos)
    assert lider_final == cod and rezagada_final == molde, "error en la replicacion"

    print("\nMoleculas hijas (replicacion semiconservativa):")
    print("  Hija 1 = hebra parental (molde) + hebra LIDER nueva")
    fila("Nueva (lider)", "5'", lider_final, "3'")
    barras()
    fila("Parental", "3'", molde, "5'")
    print("\n  Hija 2 = hebra parental (codificante) + hebra REZAGADA nueva")
    fila("Parental", "5'", cod, "3'")
    barras()
    fila("Nueva (rezagada)", "3'", rezagada_final, "5'")

    print("\nReplicacion semiconservativa completada.\n")
    return lider_final, cod