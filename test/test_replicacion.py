import os
import sys
import unittest

from Bio.Seq import Seq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import replicacion

EJEMPLO = "TTGACATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGCTTAACG"
MOLDE = str(Seq(EJEMPLO).reverse_complement())[::-1]


class TestReplicacion(unittest.TestCase):

    def test_complementaria_calculada_a_mano(self):
        self.assertEqual(replicacion.complementaria("ATGC"), "TACG")

    def test_el_cebador_es_de_arn(self):
        self.assertEqual(replicacion.cebador_arn("TTGAC"), "uugac")

    def test_las_dos_hijas_son_identicas_a_la_original(self):
        hija_1, hija_2 = replicacion.replicar_adn(EJEMPLO)
        self.assertEqual(hija_1, (EJEMPLO, MOLDE))
        self.assertEqual(hija_2, (EJEMPLO, MOLDE))

    def test_el_molde_leido_al_reves_es_el_complemento_inverso(self):
        self.assertEqual(replicacion.complementaria(EJEMPLO), MOLDE)

    def test_numero_de_fragmentos_de_okazaki(self):
        fragmentos = replicacion.sintetizar_rezagada(EJEMPLO)
        self.assertEqual(len(fragmentos), 4)
        self.assertEqual(fragmentos[-1][1], len(EJEMPLO))


if __name__ == "__main__":
    unittest.main()