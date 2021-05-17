from ufl import Argument, Coefficient
from ufl.algorithms import attach_form_arguments, strip_form_arguments


class DataCarryingArgument(Argument):
    """TODO"""
    def __init__(self, data):
        Argument.__init__(self)
        self.data = data


class DataCarryingCoefficient(Coefficient):
    """TODO"""
    def __init__(self, data):
        Coefficient.__init__(self)
        self.data = data


def test_strip_form_arguments_strips_arguments():
    data = object()
    arg = DataCarryingArgument(data)
    form = arg * dx

    stripped_form, arg_map = strip_form_arguments(form)


    # Check that the form holds a reference to the data structure
    # Check that the reference is removed

    # Check that UFL stuff still works
    assert form == stripped_form

def test_strip_form_arguments_strips_coefficients():


def test_attach_form_arguments():
    ...
