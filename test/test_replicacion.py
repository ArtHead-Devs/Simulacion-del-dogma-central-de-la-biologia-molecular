import os
import sys
import unittest

from Bio.Seq import Seq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # para encontrar replicacion.py
import replicacion

EJEMPLO = "TTGACATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGCTTAACG"


class TestReplicacion(unittest.TestCase):

    def test_complementaria_calculada_a_mano(self):
        # A-T, T-A, G-C, C-G
        self.assertEqual(replicacion.complementaria("ATGC"), "TACG")

    def test_el_cebador_es_de_arn(self):
        self.assertEqual(replicacion.cebador_arn("TTGAC"), "uugac")

    def test_las_dos_hijas_son_identicas_a_la_original(self):
        self.assertEqual(replicacion.replicar_adn(EJEMPLO), (EJEMPLO, EJEMPLO))

    def test_el_molde_leido_al_reves_es_el_complemento_inverso_de_biopython(self):
        molde = replicacion.complementaria(EJEMPLO)
        self.assertEqual(molde[::-1], str(Seq(EJEMPLO).reverse_complement()))


if __name__ == "__main__":
    unittest.main()