#!/usr/bin/env python

__authors__ = "Martin Sandve Alnes"
__date__ = "2008-03-12 -- 2008-12-02"

# Modified by Anders Logg, 2008

from ufltestcase import UflTestCase, main

from ufl import *
from ufl.algorithms import * 

# TODO: these tests only verify that the syntax is possible, how much more can we test without a form compiler?

# TODO: add more forms, covering all UFL operators

class FormsTestCase(UflTestCase):

    def test_manually_constructed_measures(self):
        # Since we can't write 'dx = dx[data]' in a non-global scope,
        # because of corner cases in the python scope rules,
        # it may be convenient to construct measures directly:
        domain_data = ('Stokes', 'Darcy')
        dx = Measure('dx')[domain_data]
        ds = Measure('ds')[domain_data]
        dS = Measure('dS')[domain_data]

        # Possible PyDOLFIN syntax:
        #ds = boundaries.dx(3) # return Measure('dx')[self](3)

    def test_measures_with_domain_data(self):
        # Configure measure with some arbitrary data object as domain_data
        domain_data = ('Stokes', 'Darcy')
        dX = dx[domain_data]

        # Build form with this domain_data
        element = FiniteElement("Lagrange", triangle, 1)
        f = Coefficient(element)
        a = f*dX(0) + f**2*dX(1)

        # Check that we get an error when using dX without domain id
        self.assertRaises(UFLException, lambda: f*dX)

        # Check that we get the right domain_data from the preprocessed form data
        fd = a.compute_form_data()
        self.assertIs(fd.domain_data['cell'], domain_data)
        self.assertIs(fd.cell_domain_data, domain_data)
        self.assertIsNone(fd.exterior_facet_domain_data)

        # Check that integral_data list is consistent as well
        f2 = f.reconstruct(count=0)
        for itd in fd.integral_data:
            t, i, itg, md = itd
            self.assertEqual(t, 'cell')
            self.assertEqual(md, {})
            self.assertEqual(itg[0].integrand(), f2**(i+1))
            self.assertIs(itg[0].measure().domain_data(), domain_data)

    def test_separated_dx(self):
        "Tests automatic summation of integrands over same domain."
        element = FiniteElement("Lagrange", triangle, 1)
        v = TestFunction(element)
        f = Coefficient(element)
        a = f*v*dx + 2*v*ds + 3*v*dx + 7*v*ds + 3*v*dx(2) + 7*v*dx(2)
        b = (f*v + 3*v)*dx + (2*v + 7*v)*ds + (3*v + 7*v)*dx(2)
        self.assertEqual(repr(a), repr(b))

    def test_source1(self):
        element = FiniteElement("Lagrange", triangle, 1)
        v = TestFunction(element)
        f = Coefficient(element)
        a = f*v*dx

    def test_source2(self):
        element = VectorElement("Lagrange", triangle, 1)
        v = TestFunction(element)
        f = Coefficient(element)
        a = dot(f,v)*dx

    def test_source3(self):
        element = TensorElement("Lagrange", triangle, 1)
        v = TestFunction(element)
        f = Coefficient(element)
        a = inner(f,v)*dx

    def test_source4(self):
        element = FiniteElement("Lagrange", triangle, 1)
        v = TestFunction(element)
        x = triangle.x
        f = sin(x[0])
        a = f*v*dx


    def test_mass1(self):
        element = FiniteElement("Lagrange", triangle, 1)
        v = TestFunction(element)
        u = TrialFunction(element)
        a = u*v*dx

    def test_mass2(self):
        element = FiniteElement("Lagrange", triangle, 1)
        v = TestFunction(element)
        u = TrialFunction(element)
        a = u*v*dx

    def test_mass3(self):
        element = VectorElement("Lagrange", triangle, 1)
        v = TestFunction(element)
        u = TrialFunction(element)
        a = dot(u,v)*dx

    def test_mass4(self):
        element = TensorElement("Lagrange", triangle, 1)
        v = TestFunction(element)
        u = TrialFunction(element)
        a = inner(u,v)*dx


    def test_stiffness1(self):
        element = FiniteElement("Lagrange", triangle, 1)
        v = TestFunction(element)
        u = TrialFunction(element)
        a = dot(grad(u), grad(v)) * dx

    def test_stiffness2(self):
        element = FiniteElement("Lagrange", triangle, 1)
        v = TestFunction(element)
        u = TrialFunction(element)
        a = inner(grad(u), grad(v)) * dx

    def test_stiffness3(self):
        element = VectorElement("Lagrange", triangle, 1)
        v = TestFunction(element)
        u = TrialFunction(element)
        a = inner(grad(u), grad(v)) * dx


    def test_stiffness_with_conductivity(self):
        velement = VectorElement("Lagrange", triangle, 1)
        telement = TensorElement("Lagrange", triangle, 1)
        v = TestFunction(velement)
        u = TrialFunction(velement)
        M = Coefficient(telement)
        a = inner(M*grad(u), grad(v)) * dx


    def test_navier_stokes(self):
        polygon = triangle
        velement = VectorElement("Lagrange", polygon, 2)
        pelement = FiniteElement("Lagrange", polygon, 1)
        TH = velement * pelement

        v, q = TestFunctions(TH)
        u, p = TrialFunctions(TH)

        f = Coefficient(velement)
        w = Coefficient(velement)
        Re = Constant(polygon)
        dt = Constant(polygon)

        a = dot(u, v) + dt*dot(dot(w, grad(u)), v)\
            - dt*Re*inner(grad(u), grad(v))\
            + dt*dot(grad(p), v)
        L = dot(f, v)
        b = dot(u, grad(q))

if __name__ == "__main__":
    main()
