import pytest
import ngsolve

def test_ParameterCF():
    p = ngsolve.fem.Parameter(23)
    assert type(p) == ngsolve.fem.Parameter
                    
