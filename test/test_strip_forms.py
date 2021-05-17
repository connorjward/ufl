"""TODO"""
import sys

from ufl import *
from ufl.algorithms import attach_form_arguments, strip_form_arguments
from ufl.core.ufl_id import attach_ufl_id
from ufl.core.ufl_type import attach_operators_from_hash_data


MIN_REF_COUNT = 2
"""The minimum value returned by sys.getrefcount. Note that this will be 3 when
run interactively."""


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

# This function is needed to provide variable scope and ensure that
# intermediate objects like the mesh and function space do not influence
# the refcount.
def _make_form(mesh_data, fs_data):
    cell = triangle
    domain = AugmentedMesh(cell, data=mesh_data)
    element = FiniteElement("Lagrange", cell, 1)
    V = AugmentedFunctionSpace(domain, element, data=fs_data)

    v = TestFunction(V) 
    u = TrialFunction(V)
    return inner(grad(v), grad(u)) * dx


def test_strip_form_arguments_strips_arguments():
    mesh_data = object()
    fs_data = object()

    # Sanity check
    assert sys.getrefcount(mesh_data) == MIN_REF_COUNT
    assert sys.getrefcount(fs_data) == MIN_REF_COUNT

    a = _make_form(mesh_data, fs_data)

    # There should now be one additional reference to both the data structures.
    assert sys.getrefcount(mesh_data) == MIN_REF_COUNT + 1
    assert sys.getrefcount(fs_data) == MIN_REF_COUNT + 1

    stripped_a, form_arg_map = strip_form_arguments(a)

    # Check that the form still works
    assert stripped_a == a

    # Delete the original form. This means that the only references to the data
    # structures are held by the mapping. This should not change the refcount.
    del a
    assert sys.getrefcount(mesh_data) == MIN_REF_COUNT + 1
    assert sys.getrefcount(fs_data) == MIN_REF_COUNT + 1

    import pdb; pdb.set_trace()

    # Check that the data references are stored in the mapping not the stripped
    # form. To do this we delete the mapping and verify that
    # the reference count for the data structures has decremented.
    del form_arg_map
    assert sys.getrefcount(mesh_data) == MIN_REF_COUNT
    assert sys.getrefcount(fs_data) == MIN_REF_COUNT


def test_strip_form_arguments_strips_coefficients():
    ...

#def # test attach
    # reattached_form = attach_form_arguments(form, mapping)

    # assert sys.getrefcount(data)

    # assert reattached_form == form

