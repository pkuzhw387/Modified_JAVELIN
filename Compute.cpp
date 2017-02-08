#include <Python.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <math.h>
#include <stdio.h>

#define WS_SIZE 1000000
#define S_SIZE 10
#define INT_KEY 4
#define CALLS 1000000
#define PI 3.1415926

/*
This is a Python-C Interface written by Haowen Zhang at Peking Univ.
The purpose is to use GSL to compute some numerical parametric integrals
in the photometric reverberation mapping.
For one-dimensional integral, gsl_integration.h is used. But for two-D ones since GSL doesn't have such functions,
I used gsl_monte.h instead. Whether the speed, accuracy and precision can meet my requirements remains to be seen.
*/

// Definition of Circle integrand. order of parameter: r, c
// r: double, radius of the disk.
// c: light speed.

static double Circle_Func(double x, void* params)
{
    double *p = static_cast<double*>(params);
    double r = p[0];
    return sqrt(r * r - x * x);
    //return r;
}



//Definition of Elliptical integrand. order of parameter: tau, i, c.
//tau: double, the value of the lag.
//i: double, the angle between the axis of the disk and the line of sight.
//c: double, light speed.

static double Ellip_Func(double x, void* params)
{
    double *p = static_cast<double*>(params);
    double tau = p[0]; double i = p[1]; double c = p[2];
    //return tau * i * c;
    return sqrt(c * c * tau * tau - x * x * (cos(i) * cos(i) + 2 * c * tau * sin(i) * x ));
}

//Definition of line-continuum integrand. Order of parameter: ti, tj, i.
//ti: double, The time on which a point in the i-th light curve was observed.
//tj: double, The time on which a point in the j-th light curve was observed.
//i: double, the angle.
//c: double, light speed.

static double LC_Func(double tau, void* params)
{
    double *p = static_cast<double*>(params);
    double ti = p[0]; double tj = p[1]; double i = p[2]; double c = p[3];
    return 2 * c * c * PI / pow(cos(i), 3) * tau * exp(-fabs(ti - tau - tj));
    //return ti*tj*i*c;
}

//Definition of line-line integrand. Order of parameters: ti, tj, i, c.
//ti: double, The time on which a point in the i-th light curve was observed.
//tj: double, The time on which a point in the j-th light curve was observed.
//i: double, the angle.
//c: double, light speed.

static double LL_Func(double x[], void* params)
{
    double *p = static_cast<double*>(params);
    double ti = p[0]; double tj = p[1]; double i = p[2]; double c = p[3];
    double tau1 = x[0]; double tau2 = x[1];
    return pow(2 * c * c * PI / pow(cos(i), 3), 2) * tau1 * tau2 * exp(-fabs(ti - tau1 - tj + tau2));

}


/*
int main()
{
    double p[4] = {3.0, 4.0, 0.4, 1.0};
    printf("%f", LC_Func(4.0, p));
    return 0;
}
*/

/* Definitions of methods */



static PyObject *
Compute_Overlap_Integrate(PyObject *self, PyObject *args)
{
    double tau, i, r, c;
    double value, value1, value2, err;

    if (!PyArg_ParseTuple(args, "ffff", &tau, &i, &r, &c))
        return NULL;

    double p1[2] = {r, c};
    double p2[3] = {tau, i, c};

    struct gsl_function_struct Circle_Integrand =
    {
        .function = Circle_Func, .params = p1
    };

    struct gsl_function_struct Ellip_Integrand =
    {
        .function = Ellip_Func, .params = p2
    };

    gsl_integration_workspace *ws1 = gsl_integration_workspace_alloc(WS_SIZE);
    gsl_integration_workspace *ws2 = gsl_integration_workspace_alloc(WS_SIZE);
    gsl_integration_qag(&Circle_Integrand, (r - c * tau) / sin(i), r, 1.0e-2, 1.0e-5,
                                WS_SIZE, INT_KEY, ws1, &value1, &err);
    gsl_integration_qag(&Ellip_Integrand, c * tau * (sin(i) - 1) / (cos(i) * cos(i)), (r - c * tau) / sin(i), 1.0e-2, 1.0e-5,
                               WS_SIZE, INT_KEY, ws2, &value2, &err);
    gsl_integration_workspace_free(ws1);
    gsl_integration_workspace_free(ws2);

    value = 2 * (value1 + value2);
    return Py_BuildValue("f", value);
}

static PyObject *
Compute_LL_Integrate_DRW(PyObject *self, PyObject *args)
{
    double ti, tj, i, c, tau0;
    double value, err;

    if(!PyArg_ParseTuple(args, "ddddd", &ti, &tj, &i, &c, &tau0))
        return NULL;
    double p[4] = {ti, tj, i, c};
    struct gsl_monte_function_struct LL_Integrand =
    {
        .f = LL_Func, .dim = 2, .params = p
    };
    double xl[2] = {0, 0};
    double xu[2] = {tau0, tau0};

    const gsl_rng_type *T;
    gsl_rng *r;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_plain_state *s = gsl_monte_plain_alloc(S_SIZE);
    gsl_monte_plain_integrate(&LL_Func, xl, xu, WS_SIZE, CALLS, r, s, &value, &err);
    gsl_monte_plain_free(s);
    gsl_rng_free(r);

    return Py_BuildValue("f", value);
}


static PyObject *
Compute_LC_Integrate_DRW(PyObject *self, PyObject *args)
{
    double ti, tj, i, c, tau0;
    double value, err;

    if(!PyArg_ParseTuple(args, "ddddd", &ti, &tj, &i, &c, &tau0))
        return NULL;

    double p[4] = {ti, tj, i, c};
    struct gsl_function_struct LC_Integrand =
    {
        .function = LC_Func, .params = p
    };

    gsl_integration_workspace *ws = gsl_integration_workspace_alloc(WS_SIZE);
    gsl_integration_qag(&LC_Integrand, 0.0, tau0, 1e-1, 1e-2, WS_SIZE, INT_KEY, ws, &value, &err);
    gsl_integration_workspace_free(ws);

    return Py_BuildValue("f", value);

}

/* Definitions of method table */

static PyMethodDef RM_Integrate_Methods[] = {

    {
        "Overlap_Integrate",
        Compute_Overlap_Integrate,
        METH_VARARGS,
        "Compute the overlap area of the disk and the isochrone."
    },

    {
        "LC_Integrate_DRW",
        Compute_LC_Integrate_DRW,
        METH_VARARGS,
        "Compute the line-continuum integrate under DRW model."
    },

    {
        "LL_Integrate_DRW",
        Compute_LL_Integrate_DRW,
        METH_VARARGS,
        "Compute the line-line integrate under DRW model."
    },

    {NULL, NULL, 0, NULL}        /* Sentinel */
};

/* Initializing function */

PyMODINIT_FUNC
initCompute(void)
{
    PyObject *m;

    m = Py_InitModule("Compute", RM_Integrate_Methods);
    if (m == NULL)
        return;
}


