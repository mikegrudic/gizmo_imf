/*
Code block containing a Brent 1973 rootfind routine.

Required initializations outside this block:
    ROOTFIND_FUNCTION(x) - #define'd macro that evaluates a function of a single
variable whose root we wish to find
    ROOTFIND_X_a, ROOTFIND_X_b - Arguments to ROOTFIND_FUNCTION that are known
to bracket the root.
    ROOTFIND_REL_X_tol - Tolerance for desired *relative* error in the root
*/

if (ROOTFUNC_a * ROOTFUNC_b > 0)
{
    PRINT_WARNING("ERROR: Bounds supplied to bracketed_roofind.h block do not bracket the root. Expanding region...");
    int fac = 1.1, iter = 0;
    while (ROOTFUNC_a * ROOTFUNC_b > 0 && iter < MAXITER)
    {
	double tmp = ROOTFIND_X_a;
        ROOTFIND_X_a = DMIN(ROOTFIND_X_a, ROOTFIND_X_b) / fac;
        ROOTFIND_X_b = DMAX(tmp, ROOTFIND_X_b) * fac;
        ROOTFUNC_a = ROOTFIND_FUNCTION(ROOTFIND_X_a);
        ROOTFUNC_b = ROOTFIND_FUNCTION(ROOTFIND_X_b);
        fac *= fac;
        iter++;
    }
    if (iter == MAXITER)
    {
        PRINT_WARNING("ERROR: Could not bracket root.\n");
    }
}

if (fabs(ROOTFUNC_a) < fabs(ROOTFUNC_b))
{ // in our convention 'a' represents
  // the bracket with the larger
  // residual
    double tmp = ROOTFUNC_a;
    ROOTFUNC_a = ROOTFUNC_b;
    ROOTFUNC_b = tmp;
    tmp = ROOTFIND_X_a;
    ROOTFIND_X_a = ROOTFIND_X_b;
    ROOTFIND_X_b = tmp;
}

double ROOTFIND_X_c = ROOTFIND_X_a, ROOTFUNC_c = ROOTFUNC_a;
int USED_BISECTION = 1, DO_BISECTION = 0, ROOTFIND_ITER = 0;
double ROOTFIND_X_c_old = ROOTFIND_X_c, ROOTFIND_X_new,
       ROOTFUNC_new = ROOTFUNC_c, ROOTFIND_REL_X_error = 1e100, EPS_TOL = 1e-8;

/* now we do a Brent 1973 method root-find */
while (ROOTFIND_REL_X_error > ROOTFIND_REL_X_tol)
{
    ROOTFIND_X_new = 0;
    if ((ROOTFUNC_a != ROOTFUNC_c) &&
        (ROOTFUNC_b != ROOTFUNC_c))
    { // inverse quadratic interpolation
        ROOTFIND_X_new += ROOTFIND_X_a * ROOTFUNC_c * ROOTFUNC_b /
                          (ROOTFUNC_a - ROOTFUNC_b) / (ROOTFUNC_a - ROOTFUNC_c);
        ROOTFIND_X_new += ROOTFIND_X_b * ROOTFUNC_c * ROOTFUNC_a /
                          (ROOTFUNC_b - ROOTFUNC_a) / (ROOTFUNC_b - ROOTFUNC_c);
        ROOTFIND_X_new += ROOTFIND_X_c * ROOTFUNC_a * ROOTFUNC_b /
                          (ROOTFUNC_c - ROOTFUNC_a) / (ROOTFUNC_c - ROOTFUNC_b);
    }
    else
    { // secant method
        ROOTFIND_X_new =
            (ROOTFIND_X_a * ROOTFUNC_b - ROOTFIND_X_b * ROOTFUNC_a) /
            (ROOTFUNC_b - ROOTFUNC_a);
    }
    DO_BISECTION = 0;
    double ROOTFIND_X_midpoint_a = 0.25 * (3 * ROOTFIND_X_a + ROOTFIND_X_b);
    if ((ROOTFIND_X_new < DMIN(ROOTFIND_X_midpoint_a, ROOTFIND_X_b)) ||
        (ROOTFIND_X_new > DMAX(ROOTFIND_X_midpoint_a, ROOTFIND_X_b)))
    {
        DO_BISECTION = 1;
    }
    if (USED_BISECTION)
    {
        if (fabs(ROOTFIND_X_new - ROOTFIND_X_b) >=
            0.5 * fabs(ROOTFIND_X_c - ROOTFIND_X_b))
        {
            DO_BISECTION = 1;
        }
        if (ROOTFIND_X_b != ROOTFIND_X_c)
        {
            if (fabs(ROOTFIND_X_b - ROOTFIND_X_c) <
                EPS_TOL * (ROOTFIND_X_b + ROOTFIND_X_c))
            {
                DO_BISECTION = 1;
            }
        }
    }
    else
    {
        if (fabs(ROOTFIND_X_new - ROOTFIND_X_b) >=
            0.5 * fabs(ROOTFIND_X_c_old - ROOTFIND_X_c))
        {
            DO_BISECTION = 1;
        }
        if (ROOTFIND_X_c_old != ROOTFIND_X_c)
        {
            if (fabs(ROOTFIND_X_c_old - ROOTFIND_X_c) <
                EPS_TOL * (ROOTFIND_X_c_old + ROOTFIND_X_c))
            {
                DO_BISECTION = 1;
            }
        }
    }
    if (DO_BISECTION)
    {
        // bisection in log space can help convergence for typical use cases; do
        // this if possible
        if ((ROOTFIND_X_b > 0) && (ROOTFIND_X_a > 0))
        {
            ROOTFIND_X_new = sqrt(ROOTFIND_X_b * ROOTFIND_X_a);
        }
        else
        {
            ROOTFIND_X_new = 0.5 * (ROOTFIND_X_b + ROOTFIND_X_a);
        }
        USED_BISECTION = 1;
    } // bisection
    else
    {
        USED_BISECTION = 0;
    }
    ROOTFUNC_new = ROOTFIND_FUNCTION(ROOTFIND_X_new);
    if (ROOTFUNC_new == 0)
    {
        break;
    }

    ROOTFIND_X_c_old = ROOTFIND_X_c;
    ROOTFIND_X_c = ROOTFIND_X_b;
    ROOTFUNC_c = ROOTFUNC_b;
    if (ROOTFUNC_a * ROOTFUNC_new < 0)
    {
        ROOTFIND_X_b = ROOTFIND_X_new;
        ROOTFUNC_b = ROOTFUNC_new;
    }
    else
    {
        ROOTFIND_X_a = ROOTFIND_X_new;
        ROOTFUNC_a = ROOTFUNC_new;
    }

    if (fabs(ROOTFUNC_a) < fabs(ROOTFUNC_b))
    {
        double tmp = ROOTFUNC_a;
        ROOTFUNC_a = ROOTFUNC_b;
        ROOTFUNC_b = tmp;
        tmp = ROOTFIND_X_a;
        ROOTFIND_X_a = ROOTFIND_X_b;
        ROOTFIND_X_b = tmp;
    }
    ROOTFIND_REL_X_error = fabs((ROOTFIND_X_b - ROOTFIND_X_a) / ROOTFIND_X_new);
    ROOTFIND_ITER++;
    if (ROOTFIND_ITER > MAXITER)
    {
        break;
    }
}

#undef ROOTFIND_FUNCTION
