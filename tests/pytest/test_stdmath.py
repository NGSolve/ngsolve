import math

import pytest
from pytest import approx
from ngsolve import *
from ngsolve import pow as ngspow
from ngsolve.fem import BSpline2D
from meshes import *


FUNCTIONS = [
    pytest.param(sqrt, 1.7,
                 lambda x: 1/(2*math.sqrt(x)),
                 lambda x: -1/(4*x**1.5), id="sqrt"),
    pytest.param(sin, 0.4, math.cos, lambda x: -math.sin(x), id="sin"),
    pytest.param(cos, 0.4, lambda x: -math.sin(x),
                 lambda x: -math.cos(x), id="cos"),
    pytest.param(tan, 0.4, lambda x: 1/math.cos(x)**2,
                 lambda x: 2*math.tan(x)/math.cos(x)**2, id="tan"),
    pytest.param(asin, 0.4, lambda x: 1/math.sqrt(1-x*x),
                 lambda x: x/(1-x*x)**1.5, id="asin"),
    pytest.param(acos, 0.4, lambda x: -1/math.sqrt(1-x*x),
                 lambda x: -x/(1-x*x)**1.5, id="acos"),
    pytest.param(atan, 0.4, lambda x: 1/(1+x*x),
                 lambda x: -2*x/(1+x*x)**2, id="atan"),
    pytest.param(sinh, 0.4, math.cosh, math.sinh, id="sinh"),
    pytest.param(cosh, 0.4, math.sinh, math.cosh, id="cosh"),
    pytest.param(exp, 0.4, math.exp, math.exp, id="exp"),
    pytest.param(log, 1.7, lambda x: 1/x, lambda x: -1/x**2, id="log"),
    pytest.param(erf, 0.4,
                 lambda x: 2/math.sqrt(math.pi)*math.exp(-x*x),
                 lambda x: -4*x/math.sqrt(math.pi)*math.exp(-x*x), id="erf"),
]


@pytest.mark.parametrize("fun,value,deriv,deriv2", FUNCTIONS)
def test_scalar_math_diff(unit_mesh_2d, fun, value, deriv, deriv2):
    ip = unit_mesh_2d(0.2, 0.2)
    u = CF(value).MakeVariable()
    f = fun(u)

    assert f.Diff(u, CF(1))(ip) == approx(deriv(value))
    assert f.Diff(u)(ip) == approx(deriv(value))
    assert f.Diff(u, CF(1)).Diff(u, CF(1))(ip) == approx(deriv2(value))


@pytest.mark.parametrize("fun,value,deriv,deriv2", FUNCTIONS)
def test_vector_math_diff(unit_mesh_2d, fun, value, deriv, deriv2):
    ip = unit_mesh_2d(0.2, 0.2)
    values = (value, value+0.2)
    direction = (2, -3)
    u = CF(values).MakeVariable()
    f = fun(u)

    expected = tuple(deriv(x)*d for x, d in zip(values, direction))
    df = f.Diff(u, CF(direction))
    assert tuple(df.dims) == (2,)
    assert df(ip) == approx(expected)
    jac2 = df.Diff(u)
    assert tuple(jac2.dims) == (2, 2)
    assert jac2(ip) == approx((deriv2(values[0])*direction[0], 0,
                               0, deriv2(values[1])*direction[1]))
    assert df.Diff(df)(ip) == approx((1, 0, 0, 1))

    jac = f.Diff(u)
    assert tuple(jac.dims) == (2, 2)
    assert jac(ip) == approx((deriv(values[0]), 0, 0, deriv(values[1])))


@pytest.mark.parametrize("dims", [(2,), (2, 2)])
def test_componentwise_diff_replace(unit_mesh_2d, dims):
    ip = unit_mesh_2d(0.2, 0.2)
    n = math.prod(dims)
    values = (0.2, 0.4, 0.6, 0.8)[:n]
    replacements = (0.3, 0.5, 0.7, 0.9)[:n]
    direction = (2, 3, 5, 7)[:n]
    u = CF(values, dims=dims).MakeVariable()
    v = CF(replacements, dims=dims).MakeVariable()
    df = sin(u).Diff(u, CF(direction, dims=dims))
    replaced = df.Replace({u: v})
    expected = tuple(math.cos(x)*d for x, d in zip(replacements, direction))

    assert tuple(replaced.dims) == dims
    assert replaced(ip) == approx(expected)
    assert replaced.Compile()(ip) == approx(expected)


def test_atan2_diff(unit_mesh_2d):
    ip = unit_mesh_2d(0.2, 0.2)
    yval, xval = 0.7, 1.2
    y = CF(yval).MakeVariable()
    x = CF(xval).MakeVariable()
    f = atan2(y, x)
    r2 = xval*xval + yval*yval

    assert f.Diff(y)(ip) == approx(xval/r2)
    assert f.Diff(x)(ip) == approx(-yval/r2)
    assert f.Diff(y, CF(1)).Diff(y, CF(1))(ip) == approx(-2*xval*yval/r2**2)
    assert f.Diff(y, CF(1)).Diff(x, CF(1))(ip) == approx((yval*yval-xval*xval)/r2**2)


def test_atan2_vector_diff(unit_mesh_2d):
    ip = unit_mesh_2d(0.2, 0.2)
    yval, xval, direction = (0.7, 0.4), (1.2, 0.9), (2, -3)
    y = CF(yval).MakeVariable()
    x = CF(xval).MakeVariable()
    f = atan2(y, x)
    r2 = tuple(a*a+b*b for a, b in zip(yval, xval))

    expected = tuple(d*b/r for d, b, r in zip(direction, xval, r2))
    df = f.Diff(y, CF(direction))
    assert tuple(df.dims) == (2,)
    assert df(ip) == approx(expected)

    jac = f.Diff(y)
    assert tuple(jac.dims) == (2, 2)
    assert jac(ip) == approx((xval[0]/r2[0], 0, 0, xval[1]/r2[1]))

    expected = tuple(-d*a/r for d, a, r in zip(direction, yval, r2))
    assert f.Diff(x, CF(direction))(ip) == approx(expected)
    jac = f.Diff(x)
    assert tuple(jac.dims) == (2, 2)
    assert jac(ip) == approx((-yval[0]/r2[0], 0, 0, -yval[1]/r2[1]))


@pytest.mark.parametrize("value,exponent,expected", [(-2.0, 3.0, 12.0),
                                                       (0.0, 1.0, 1.0)])
def test_pow_constant_exponent_diff(unit_mesh_2d, value, exponent, expected):
    ip = unit_mesh_2d(0.2, 0.2)
    u = CF(1).MakeVariable()
    x = u+(value-1)
    f = ngspow(x, exponent)

    assert f.Diff(u, CF(1))(ip) == approx(expected)
    assert f.Diff(u)(ip) == approx(expected)


def test_pow_vector_diff(unit_mesh_2d):
    ip = unit_mesh_2d(0.2, 0.2)
    baseval, exponentval, direction = (2.0, 3.0), (0.4, 1.2), (2, -3)
    base = CF(baseval).MakeVariable()
    exponent = CF(exponentval).MakeVariable()
    f = ngspow(base, exponent)

    expected = tuple(d*b*a**(b-1)
                     for a, b, d in zip(baseval, exponentval, direction))
    df = f.Diff(base, CF(direction))
    assert tuple(df.dims) == (2,)
    assert df(ip) == approx(expected)

    jac = f.Diff(base)
    assert tuple(jac.dims) == (2, 2)
    assert jac(ip) == approx((exponentval[0]*baseval[0]**(exponentval[0]-1), 0,
                              0, exponentval[1]*baseval[1]**(exponentval[1]-1)))

    expected = tuple(d*a**b*math.log(a)
                     for a, b, d in zip(baseval, exponentval, direction))
    assert f.Diff(exponent, CF(direction))(ip) == approx(expected)
    jac = f.Diff(exponent)
    assert tuple(jac.dims) == (2, 2)
    assert jac(ip) == approx((baseval[0]**exponentval[0]*math.log(baseval[0]), 0,
                              0, baseval[1]**exponentval[1]*math.log(baseval[1])))


def test_bspline_diff(unit_mesh_2d):
    ip = unit_mesh_2d(0.2, 0.2)
    spline = BSpline(3, [0, 1, 2, 3, 4], [0, 1, 4, 9, 16])
    dspline = spline.Differentiate()
    values, direction = (2.3, 2.7), (2, -3)
    u = CF(values).MakeVariable()
    f = spline(u)

    expected = tuple(dspline(x)*d for x, d in zip(values, direction))
    df = f.Diff(u, CF(direction))
    assert tuple(df.dims) == (2,)
    assert df(ip) == approx(expected)

    jac = f.Diff(u)
    assert tuple(jac.dims) == (2, 2)
    assert jac(ip) == approx((dspline(values[0]), 0, 0, dspline(values[1])))


def test_bspline2d_diff(unit_mesh_2d):
    spline = BSpline2D([0, 1], [0, 1], [0, 3, 2, 5])
    x = CF(0.2).MakeVariable()
    y = CF(0.4).MakeVariable()
    f = spline(x, y)

    for var in [x, y]:
        with pytest.raises(Exception, match="does not provide a derivative"):
            f.Diff(var)
