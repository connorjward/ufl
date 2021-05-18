from ufl import *
from ufl.algorithms import attach_form_arguments, strip_form_arguments
from ufl.core.ufl_id import attach_ufl_id
from ufl.core.ufl_type import attach_operators_from_hash_data


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


def test_strip_form_arguments():
    mesh_data = "MESH DATA"
    fs_data = "FUNCTION SPACE DATA"

    cell = triangle
    domain = AugmentedMesh(cell, data=mesh_data)
    element = FiniteElement("Lagrange", cell, 1)
    V = AugmentedFunctionSpace(domain, element, data=fs_data)

    v = TestFunction(V) 
    u = TrialFunction(V)
    f = Coefficient(V)
    k = Constant(V)
    form = k * f * inner(grad(v), grad(u)) * dx

    stripped_form, mapping = strip_form_arguments(form)

    # Check that the form still works
    # TODO: How can I get this to work?
    # assert stripped_form == form

    # Verify the contents of the mapping
    for old_const, new_const in zip(form.constants(), stripped_form.constants()):
        assert not hasattr(new_const.ufl_domain(), "data")
        assert mapping[new_const] is old_const

    for old_arg, new_arg in zip(form.arguments(), stripped_form.arguments()):
        assert not hasattr(new_arg.ufl_function_space(), "data")
        assert not hasattr(new_arg.ufl_function_space().ufl_domain(), "data")
        assert mapping[new_arg] is old_arg

    for old_coeff, new_coeff in zip(form.coefficients(), stripped_form.coefficients()):
        assert not hasattr(new_coeff.ufl_function_space(), "data")
        assert not hasattr(new_coeff.ufl_function_space().ufl_domain(), "data")
        assert mapping[new_coeff] is old_coeff

    reattached_form = attach_form_arguments(stripped_form, mapping)

    # Check that the form still works
    # TODO: How?
    # assert reattached_form == form

    # Check that the arguments have been replaced correctly
    for old_const, new_const in zip(form.constants(), reattached_form.constants()):
        assert old_const is new_const
    for old_arg, new_arg in zip(form.arguments(), reattached_form.arguments()):
        assert old_arg is new_arg
    for old_coeff, new_coeff in zip(form.coefficients(), reattached_form.coefficients()):
        assert old_coeff is new_coeff

