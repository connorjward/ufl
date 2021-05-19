import gc
import sys

import pytest

from ufl import *
from ufl.algorithms import attach_form_arguments, strip_form_arguments
from ufl.core.ufl_id import attach_ufl_id
from ufl.core.ufl_type import attach_operators_from_hash_data


MIN_REF_COUNT = 2
"""The minimum value returned by sys.getrefcount."""


@attach_operators_from_hash_data
@attach_ufl_id
class AugmentedMesh(Mesh):
    def __init__(self, *args, data):
        super().__init__(*args)
        self.data = data


@attach_operators_from_hash_data
class AugmentedFunctionSpace(FunctionSpace):
    def __init__(self, *args, data):
        super().__init__(*args)
        self.data = data


class AugmentedCoefficient(Coefficient):
    def __init__(self, *args, data):
        super().__init__(*args)
        self.data = data


def _make_form(mesh_data, fs_data, coeff_data):
    cell = triangle
    domain = AugmentedMesh(cell, data=mesh_data)
    element = FiniteElement("Lagrange", cell, 1)
    V = AugmentedFunctionSpace(domain, element, data=fs_data)

    v = TestFunction(V) 
    u = TrialFunction(V)
    f = AugmentedCoefficient(V, data=coeff_data)
    k = Constant(V)

    return k*f*inner(grad(v), grad(u))*dx


def test_strip_form_arguments_strips_data_refs():
    mesh_data = object()
    fs_data = object()
    coeff_data = object()

    assert sys.getrefcount(mesh_data) == MIN_REF_COUNT
    assert sys.getrefcount(fs_data) == MIN_REF_COUNT
    assert sys.getrefcount(coeff_data) == MIN_REF_COUNT

    form = _make_form(mesh_data, fs_data, coeff_data)
    stripped_form, mapping = strip_form_arguments(form)

    assert sys.getrefcount(mesh_data) == MIN_REF_COUNT + 1
    assert sys.getrefcount(fs_data) == MIN_REF_COUNT + 1
    assert sys.getrefcount(coeff_data) == MIN_REF_COUNT + 1

    del form, mapping
    gc.collect()  # This is needed to update the refcounts

    assert sys.getrefcount(mesh_data) == MIN_REF_COUNT
    assert sys.getrefcount(fs_data) == MIN_REF_COUNT
    assert sys.getrefcount(coeff_data) == MIN_REF_COUNT


def test_strip_form_arguments_does_not_change_form():
    form = _make_form(None, None, None)
    stripped_form, mapping = strip_form_arguments(form)

    assert stripped_form.signature() == form.signature()
    assert attach_form_arguments(stripped_form, mapping) == form
