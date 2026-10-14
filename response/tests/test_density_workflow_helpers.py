#!/usr/bin/env python3
"""Pure helper tests; no student or many-body solver is launched."""
from pathlib import Path
import sys
import unittest

import numpy as np

HERE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(HERE))
from mmhk_density_workflow import choose_mu,k_grid,finite_ring_probe_weights,number_minimizers


class DensityHelpers(unittest.TestCase):
    def test_supported_q_is_unit_weight(self):
        for m in range(6):
            actual=finite_ring_probe_weights(2*np.pi*m/6,6)
            expected=np.zeros(6);expected[m]=1
            np.testing.assert_allclose(actual,expected,atol=1e-14)

    def test_window_parseval(self):
        for q in [-5.,-.1,0.,.62831853,np.pi,9.1]:
            weights=finite_ring_probe_weights(q,6)
            self.assertTrue(np.all(weights>=0))
            self.assertAlmostEqual(weights.sum(),1.)

    def test_grid_conventions(self):
        student,w=k_grid(6,16,'student');midpoint,v=k_grid(6,16,'midpoint')
        self.assertAlmostEqual(student[0],-np.pi/6);self.assertAlmostEqual(student[-1],np.pi/6)
        self.assertGreater(midpoint[0],-np.pi/6);self.assertLess(midpoint[-1],np.pi/6)
        np.testing.assert_allclose(student,-student[::-1]);np.testing.assert_allclose(midpoint,-midpoint[::-1])
        self.assertAlmostEqual(w.sum(),1.);self.assertAlmostEqual(v.sum(),1.)
        for count in [4,5,16]:
            periodic,weights=k_grid(6,count,'periodic')
            self.assertEqual(len(np.unique(periodic)),count)
            self.assertTrue(np.any(periodic==0))
            self.assertLess(periodic[-1],np.pi/6)
            self.assertAlmostEqual(weights.sum(),1.)

    def test_discontinuous_density_is_not_fabricated(self):
        # L=1 toy envelope: particle sectors N=0,1,2. N=1 is unstable.
        energies=np.array([[0.,2.,1.]])
        result=choose_mu(energies,np.array([1.]),.45,1)
        self.assertEqual(result['filling'],0.)
        self.assertFalse(result['exact_target_attained'])
        self.assertEqual(number_minimizers(energies,.5),[[0,2]])

    def test_plateau_mu_matches_brute_candidates(self):
        file=HERE/'results/physical_mmhk_L6_validation/U8/ground_energy_table.npz'
        with np.load(file,allow_pickle=False) as d:energies=d['E0_by_N'];weights=d['weights']
        result=choose_mu(energies,weights,.45,6)
        crossings=[]
        for row in energies:
            for a in range(13):
                for b in range(a+1,13):crossings.append((row[b]-row[a])/(b-a))
        crossings=np.unique(np.round(crossings,8));candidates=[]
        for left,right in zip(crossings[:-1],crossings[1:]):
            if right-left<1e-7:continue
            mu=(left+right)/2;numbers=np.argmin(energies-mu*np.arange(13),axis=1)
            candidates.append(float(weights@numbers/12))
        self.assertAlmostEqual(abs(result['filling']-.45),min(abs(n-.45) for n in candidates))


if __name__=='__main__':unittest.main()
