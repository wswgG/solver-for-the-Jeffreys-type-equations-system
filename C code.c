/*** STRUCT ****************************************************************
 *              The following struct holds the equation parameters         *
 ***************************************************************************/

typedef struct EqParms {
	double *R;           /* strength of each promoter--always >= 0. */
	double *T;           /* the genetic interconnect matrix */
	double *E;           /* the external input regulatory matrix */
	double *M;          /* the maternal (constant in time) input regulatory matrix */
	double *P;          /* the couples input regulatory matrix */
	double *h;           /* reg. coeff. for generic TFs on synthesis of gene */
	double *d;           /* spatial interaction at gastrulation--always >= 0. */
	double *lambda; /*protein half lives--always >= 0. */
	double *tau;        /* delay times for the proteins */
	double *m;        /* delay times for the proteins */
	double *mm;        /* delay times for the proteins */
} EqParms;


/***    exFinDif: propagates vin (of size n) from tin to tout by the      ***
 *     explicit Finite-Difference method of 1st or 2nd order accuracy;      *
 *                       the result is returned by vout                     *
 *                                                                          *
 * ************************************************************************ *
 *                                                                          *
 * NOTATIONS :                                                              *
 * vin        <--   input meanings of concentration                         *
 * win        <--   input meanings of concentration function gradient       *
 * vout       <--   output meanings of concentration                        *
 * wout       <--   output meanings of concentration function gradient      *
 * tin        <--   start time                                              *
 * tout       <--   end time                                                *
 * stephint   <--   stepsize                                                *
 * accuracy   <--   accuracy for adaptive stepsize solvers                  *
 * n          <--   length of vin, win, vout, wout                          *
 *                                                                          *
 *** ******************************************************************* ***/

ODESolutionPoly* exFinDif(double *vin, double *win, double *vout, double *wout, double tin, double tout,
	 double stephint, double accuracy, int n, FILE *slog)
{
	double **v;               /* three time layers of solution */
  double *temp;             /* temporary variable */
	double theta;             /* time step size */
	double t;                 /* current time moment */
	double m;         		  /* used to calculate time step */
	int flag;                 /* help variable for calculating time step */
	int timeNum;              /* time layer index */
	int i;                    /* help indexes */
	int M;                    /* time layer count */
	ODESolutionPoly *Solution;
   /* the do-nothing case: too small steps dealt with under usual */
	if (tin == tout) {
		return NULL;
	}
	Solution = ODESolutionPoly_createWith(tin, tout, n);
	/* allocating memory */
	v = (double**) malloc(3 * sizeof(double*));
	v[0] = (double*) malloc(n * sizeof(double));
	v[1] = (double*) malloc(n * sizeof(double));
	v[2] = (double*) malloc(n * sizeof(double));
  /* calculating time step */
	m = floor(fabs(tout - tin) / stephint + 0.5);  /* see comment on stephint above */
	if (m < 1.) {
		m = 1.;                                /* we'll have to do at least one step */
	}
	theta = (tout - tin) / m;                 /* real stepsize calculated here */
	M = (int)m;                              /* int number of steps */
  t = tin;                                 /* set current time */
  if (t == t + theta) {
    error("Finite Difference: stephint of %g too small!", theta);
	}
  /* filling solution vectors by initial values */
	for (i = 0; i < n; i++) {
		v[0][i] = win[i];
		v[1][i] = vin[i];
		v[2][i] = 0.0;
	}
    /* main loop for time steps */
    /* for (timeNum = 0; timeNum < M; timeNum++) {
      v[2][0] = theta;
      v[2][1] = tin;
      p_deriv2(v[1], v[0], t, v[2], n);
      t += theta;
      ODESolutionPoly_addtoInterpolating(Solution, v[2], t);
      // preparing for calculating next time layer solution
      for (i = 0; i < n; i++) {
         v[0][i] = v[1][i];
         v[1][i] = v[2][i];
         v[2][i] = 0.0;
      }
   } */

   for (timeNum = 0; timeNum < M - 2; timeNum++) {
     v[2][0] = theta;
     v[2][1] = tin;
     p_deriv2(v[1], v[0], t, v[2], n);
     t += theta;
     ODESolutionPoly_addtoInterpolating(Solution, v[2], t);
     // preparing for calculating next time layer solution
     temp = v[0];
     v[0] = v[1];
     v[1] = v[2];
     v[2] = temp;
   }

   // last but one time step
   free(v[2]);
   v[2] = wout;
   v[2][0] = theta;
   v[2][1] = tin;
   p_deriv2(v[1], v[0], t, v[2], n);
   t += theta;
   ODESolutionPoly_addtoInterpolating(Solution, v[2], t);
   // preparing for calculating last time layer solution
   free(v[0]);
   v[0] = v[1];
   v[1] = v[2];
   v[2] = vout;

   // last time step
   v[2][0] = theta;
   v[2][1] = tin;
   p_deriv2(v[1], v[0], t, v[2], n);
   t += theta;
   ODESolutionPoly_addtoInterpolating(Solution, v[2], t);
   
   ODESolutionPoly_makeInterpolating(Solution);
   free(v[0]);
   free(v);
   return Solution;
}


/***                       Dvdt31stJeffDirichlet                         ***
 *                                                                         *
 *      uses explicit Finite-Difference Method of 1st order accuracy       *
 *       for the Jeffreys type equation, calculates next time step         *
 *                Dirichlet boundary conditions are used                   *
 *                                                                         *
 * *********************************************************************** *
 *                                                                         *
 * NOTATIONS :                                                             *
 * v1        <--   input meanings of concentration                         *
 * v0        <--   input meanings of concentration function gradient       *
 * vdot      <--   output meanings of concentration                        *
 * t         <--   current time moment                                     *
 * n         <--   length of v1, v0, vdot                                  *
 *                                                                         *
 *** ****************************************************************** ***/
void Dvdt31stJeffDirichlet(double *v1, double *v0, double t, double *vdot, int n)
{
  static unsigned int samecycle = 0;     /* same cleavage cycle as before? */
  double vinput1, vinput3;         /* help variables */
	double theta = vdot[0];          /* time step size */
  double tin = vdot[1];            /* initial time   */
	double g;                        /* g-function value */
	double gDer;                     /* g-function derivative value */
	double eps;                      /* help cofficient */
	double alpha;                    /* help cofficient */
	double beta;                     /* help cofficient */
	double chi;                      /* help cofficient */
  double delta;                    /* help cofficient */
	double prevInpT, curInpT;        /* inputs from considered genes */
	double prevInpE, curInpE;        /* inputs from external genes */
	double inpM;                     /* input from maternal genes */
  double g10 = vdot[3] - vdot[2];
  double g1L = vdot[5] - vdot[4];
  int i, j, k, ap;                 /* local indexes */
  
	if (ccycle != samecycle) {
		GetD(t, lparm->d, defs->diff_schedule, D);
		GetL(t, lparm->m, defs->mob_schedule, M);
		/* GetLL(t, lparm->mm, defs->mob_schedule, MM); */
    GetD(t, lparm->mm, defs->diff_schedule, MM);
		GetRCoef(t, &mit_factor);
    vmat = gcmd_get_maternal_inputs(ccycle);
		samecycle = ccycle;
	}
	if (rule == INTERPHASE) {
	  if (defs->egenes > 0) {
	    gcmd_get_external_inputs(t, vext, num_nucs * defs->egenes);
			gcdm_get_curveSlopeExtInputs(t, vext2, num_nucs * defs->egenes);
		}
		for (ap = 1; ap < num_nucs - 1; ap++) {
			for (i = 0; i < defs->ngenes; i++) {
				chi = 1.0 / (1 + lparm->lambda[i] * lparm->tau[i]);       
				eps = lparm->tau[i] * chi;
				alpha = D[i] * chi;
				beta = -lparm->lambda[i] * chi;
				delta = -lparm->tau[i] * MM[i] * chi;
				prevInpT = 0.0;
				curInpT = 0.0;
				prevInpE = 0.0;
				curInpE = 0.0;
				inpM = 0.0;
				for (k = 0; k < defs->ngenes; k++) {
					curInpT += lparm->T[defs->ngenes * i + k] * v1[defs->ngenes * ap + k];
					prevInpT += lparm->T[defs->ngenes * i + k] * v0[defs->ngenes * ap + k];
				}
				for (k = 0; k < defs->egenes; k++) {
               curInpE += lparm->E[defs->egenes * i + k] * vext[defs->egenes * ap + k];
					prevInpE += lparm->E[defs->egenes * i + k] * vext2[defs->egenes * ap + k];
				}
				for (k = 0; k < defs->mgenes; k++) {
					inpM += lparm->M[defs->mgenes * i + k] * vmat[defs->mgenes * ap + k];
				}
				vinput1 = curInpT + curInpE + inpM + lparm->h[i];
				vinput3 = 1 / sqrt(1 + vinput1 * vinput1);
				if (gofu == Tanh) {
					g = lparm->R[i] / 2.0 * (1.0 + tanh(vinput1));
					gDer = lparm->R[i] / (2.0 * cosh(vinput1) * cosh(vinput1));
				} else if (gofu == Sqrt) {
					g = 0.5 * lparm->R[i] * (1.0 + vinput1 * vinput3);
					gDer = 0.5 * lparm->R[i] * vinput3 * vinput3 * vinput3;
				} else if (gofu == Lin) {
				    if (vinput1 > -1 && vinput1 < 1) {
						g = lparm->R[i] * (vinput1 + 1.0) / 2.0;
						gDer = lparm->R[i] / 2.0;
					}
					else {
						gDer = 0.0;
						if (vinput1 > 1)
							g = lparm->R[i];
						else
							g = 0.0;
					}
				} else {
               error("Dvdt31stJeffDirichlet: unknown g(u)");
				}
        /* calculating next time step */
        vdot[defs->ngenes * ap + i] =
          v1[defs->ngenes * ap + i] * (2 * eps / (theta * theta) + 2 * delta / theta - 2 * alpha + beta) +
          eps * gDer / theta * (curInpT - prevInpT) + eps * gDer * prevInpE + chi * g +
          (alpha - delta / theta) * (v1[defs->ngenes * (ap + 1) + i] + v1[defs->ngenes * (ap - 1) + i]) +
          delta / theta * (v0[defs->ngenes * (ap + 1) + i] + v0[defs->ngenes * (ap - 1) + i]) +
          v0[defs->ngenes * ap + i] * (-eps / (theta * theta) + 1 / (2.0 * theta) - 2 * delta / theta);

          vdot[defs->ngenes * ap + i] /= eps / (theta * theta) + 1 / (2.0 * theta);
         }
      }

      /* Dirichlet boundary conditions */
      for (i = 0; i < defs->ngenes; i++) {
        vdot[i] = vext[(defs->egenes - defs->ngenes) + i];
        vdot[defs->ngenes * (num_nucs - 1) + i] = vext[defs->egenes * (num_nucs - 1) + (defs->egenes - defs->ngenes) + i];
      }
   } else if (rule == MITOSIS) {
      for (ap = 1; ap < num_nucs - 1; ap++) {
        for (i = 0; i < defs->ngenes; i++) {
          chi = 1.0 / (1 + lparm->lambda[i] * lparm->tau[i]);       
  				eps = lparm->tau[i] * chi;
  				alpha = D[i] * chi;
  				beta = -lparm->lambda[i] * chi;
  				delta = -lparm->tau[i] * MM[i] * chi;

          /* calculating next time step */
          vdot[defs->ngenes * ap + i] =
            v1[defs->ngenes * ap + i] * (2 * eps / (theta * theta) + 2 * delta / theta - 2 * alpha + beta) +
            (alpha - delta / theta) * (v1[defs->ngenes * (ap + 1) + i] + v1[defs->ngenes * (ap - 1) + i]) +
            delta / theta * (v0[defs->ngenes * (ap + 1) + i] + v0[defs->ngenes * (ap - 1) + i]) +
            v0[defs->ngenes * ap + i] * (-eps / (theta * theta) + 1 / (2.0 * theta) - 2 * delta / theta);
            
            vdot[defs->ngenes * ap + i] /= eps / (theta * theta) + 1 / (2.0 * theta);
         }
      }

      /* Dirichlet boundary conditions */
      for (i = 0; i < defs->ngenes; i++) {
        vdot[i] = vext[(defs->egenes - defs->ngenes) + i];
        vdot[defs->ngenes * (num_nucs - 1) + i] = vext[defs->egenes * (num_nucs - 1) + (defs->egenes - defs->ngenes) + i];
      }
   } else {
      error("Dvdt31stJeffDirichlet: Bad rule %i sent to Dvdt31stJeffDirichlet", rule);
	}
}

/***                       Dvdt31stJeffNeumann                           ***
 *                                                                         *
 *      uses explicit Finite-Difference Method of 1st order accuracy       *
 *       for the Jeffreys type equation, calculates next time step         *
 *                  Neumann boundary conditions are used                   *
 *                                                                         *
 * *********************************************************************** *
 *                                                                         *
 * NOTATIONS :                                                             *
 * v1        <--   input meanings of concentration                         *
 * v0        <--   input meanings of concentration function gradient       *
 * vdot      <--   output meanings of concentration                        *
 * t         <--   current time moment                                     *
 * n         <--   length of v1, v0, vdot                                  *
 *                                                                         *
 *** ****************************************************************** ***/
void Dvdt31stJeffNeumann(double *v1, double *v0, double t, double *vdot, int n)
{
  static unsigned int samecycle = 0;     /* same cleavage cycle as before? */
	double vinput1, vinput3;         /* help variables */
	double theta = vdot[0];          /* time step size */
   double tin = vdot[1];            /* initial time   */
	double g;                        /* g-function value */
	double gDer;                     /* g-function derivative value */
	double eps;                      /* help cofficient */
	double alpha;                    /* help cofficient */
	double beta;                     /* help cofficient */
	double chi;                      /* help cofficient */
   double delta;                    /* help cofficient */
	double prevInpT, curInpT;        /* inputs from considered genes */
	double prevInpE, curInpE;        /* inputs from external genes */
	double inpM;                     /* input from maternal genes */
   double g10 = vdot[3] - vdot[2];
   double g1L = vdot[5] - vdot[4];
   int i, j, k, ap;                 /* local indexes */
  
	if (ccycle != samecycle) {
		GetD(t, lparm->d, defs->diff_schedule, D);
		GetL(t, lparm->m, defs->mob_schedule, M);
		/* GetLL(t, lparm->mm, defs->mob_schedule, MM); */
      GetD(t, lparm->mm, defs->diff_schedule, MM);
		GetRCoef(t, &mit_factor);
      vmat = gcmd_get_maternal_inputs(ccycle);
		samecycle = ccycle;
	}
	if (rule == INTERPHASE) {
      if (defs->egenes > 0) {
			gcmd_get_external_inputs(t, vext, num_nucs * defs->egenes);
			gcdm_get_curveSlopeExtInputs(t, vext2, num_nucs * defs->egenes);
		}
		for (ap = 0; ap < num_nucs; ap++) {
			for (i = 0; i < defs->ngenes; i++) {
				chi = 1.0 / (1 + lparm->lambda[i] * lparm->tau[i]);       
				eps = lparm->tau[i] * chi;
				alpha = D[i] * chi;
				beta = -lparm->lambda[i] * chi;
            delta = -lparm->tau[i] * MM[i] * chi;
				prevInpT = 0.0;
				curInpT = 0.0;
				prevInpE = 0.0;
				curInpE = 0.0;
				inpM = 0.0;
				for (k = 0; k < defs->ngenes; k++) {
					curInpT += lparm->T[defs->ngenes * i + k] * v1[defs->ngenes * ap + k];
					prevInpT += lparm->T[defs->ngenes * i + k] * v0[defs->ngenes * ap + k];
				}
				for (k = 0; k < defs->egenes; k++) {
               curInpE += lparm->E[defs->egenes * i + k] * vext[defs->egenes * ap + k];
					prevInpE += lparm->E[defs->egenes * i + k] * vext2[defs->egenes * ap + k];
				}
				for (k = 0; k < defs->mgenes; k++) {
					inpM += lparm->M[defs->mgenes * i + k] * vmat[defs->mgenes * ap + k];
				}
				vinput1 = curInpT + curInpE + inpM + lparm->h[i];
				vinput3 = 1 / sqrt(1 + vinput1 * vinput1);
				if (gofu == Tanh) {
					g = lparm->R[i] / 2.0 * (1.0 + tanh(vinput1));
					gDer = lparm->R[i] / (2.0 * cosh(vinput1) * cosh(vinput1));
				} else if (gofu == Sqrt) {
					g = 0.5 * lparm->R[i] * (1.0 + vinput1 * vinput3);
					gDer = 0.5 * lparm->R[i] * vinput3 * vinput3 * vinput3;
				} else if (gofu == Lin) {
					if (vinput1 > -1 && vinput1 < 1) {
						g = lparm->R[i] * (vinput1 + 1.0) / 2.0;
						gDer = lparm->R[i] / 2.0;
					}
					else {
						gDer = 0.0;
						if (vinput1 > 1)
							g = lparm->R[i];
						else
							g = 0.0;
					}
				} else {
               error("Dvdt31stJeffNeumann: unknown g(u)");
				}

            /* calculating next time step */
            vdot[defs->ngenes * ap + i] =
               v1[defs->ngenes * ap + i] * (2 * eps / (theta * theta) + 2 * delta / theta- 2 * alpha + beta) +
               eps * gDer / theta * (curInpT - prevInpT) + eps * gDer * prevInpE + chi * g +
               (alpha - delta / theta) * (v1[defs->ngenes * (ap + 1) + i] + v1[defs->ngenes * (ap - 1) + i]) +
               delta / theta * (v0[defs->ngenes * (ap + 1) + i] + v0[defs->ngenes * (ap - 1) + i]) +
               v0[defs->ngenes * ap + i] * (-eps / (theta * theta) + 1 / (2.0 * theta) - 2 * delta / theta);
            vdot[defs->ngenes * ap + i] /= eps / (theta * theta) + 1 / (2.0 * theta);
         }

         /* Neumann boundary conditions */
         for (i = 0; i < defs->ngenes; i++) {
            // the far left nucleus
            vdot[i] =
               v1[i] * (2 * eps / (theta * theta) + 2 * delta / theta- 2 * alpha + beta) +
               eps * gDer / theta * (curInpT - prevInpT) + eps * gDer * prevInpE + chi * g +
               2 * (alpha - delta / theta) * v1[defs->ngenes + i] +
               2 * delta / theta * v0[defs->ngenes + i] +
               v0[i] * (-eps / (theta * theta) + 1 / (2.0 * theta) - 2 * delta / theta);
            vdot[i] /= eps / (theta * theta) + 1 / (2.0 * theta);

            // the far right nucleus
            vdot[defs->ngenes * (num_nucs - 1) + i] =
               v1[defs->ngenes * (num_nucs - 1) + i] * (2 * eps / (theta * theta) + 2 * delta / theta- 2 * alpha + beta) +
               eps * gDer / theta * (curInpT - prevInpT) + eps * gDer * prevInpE + chi * g +
               2 * (alpha - delta / theta) * v1[defs->ngenes * (num_nucs - 2) + i] +
               2 * delta / theta * v0[defs->ngenes * (num_nucs - 2) + i] +
               v0[defs->ngenes * (num_nucs - 1) + i] * (-eps / (theta * theta) + 1 / (2.0 * theta) - 2 * delta / theta);
            vdot[defs->ngenes * (num_nucs - 1) + i] /= eps / (theta * theta) + 1 / (2.0 * theta);
         }
      }
   } else if (rule == MITOSIS) {
      for (ap = 0; ap < num_nucs; ap++) {
         for (i = 0; i < defs->ngenes; i++) {
            chi = 1.0 / (1 + lparm->lambda[i] * lparm->tau[i]);
            eps = lparm->tau[i] * chi;
            alpha = D[i] * chi;
            beta = -lparm->lambda[i] * chi;
            delta = -lparm->tau[i] * MM[i] * chi;

            /* calculating next time step */
            vdot[defs->ngenes * ap + i] =
               v1[defs->ngenes * ap + i] * (2 * eps / (theta * theta) + 2 * delta / theta- 2 * alpha + beta) +
               (alpha - delta / theta) * (v1[defs->ngenes * (ap + 1) + i] + v1[defs->ngenes * (ap - 1) + i]) +
               delta / theta * (v0[defs->ngenes * (ap + 1) + i] + v0[defs->ngenes * (ap - 1) + i]) +
               v0[defs->ngenes * ap + i] * (-eps / (theta * theta) + 1 / (2.0 * theta) - 2 * delta / theta);
            vdot[defs->ngenes * ap + i] /= eps / (theta * theta) + 1 / (2.0 * theta);
         }

         /* Neumann boundary conditions */
         for (i = 0; i < defs->ngenes; i++) {
            // the far left nucleus
            vdot[i] =
               v1[i] * (2 * eps / (theta * theta) + 2 * delta / theta- 2 * alpha + beta) +
               2 * (alpha - delta / theta) * v1[defs->ngenes + i] +
               2 * delta / theta * v0[defs->ngenes + i] +
               v0[i] * (-eps / (theta * theta) + 1 / (2.0 * theta) - 2 * delta / theta);
            vdot[i] /= eps / (theta * theta) + 1 / (2.0 * theta);

            // the far right nucleus
            vdot[defs->ngenes * (num_nucs - 1) + i] =
               v1[defs->ngenes * (num_nucs - 1) + i] * (2 * eps / (theta * theta) + 2 * delta / theta- 2 * alpha + beta) +
               2 * (alpha - delta / theta) * v1[defs->ngenes * (num_nucs - 2) + i] +
               2 * delta / theta * v0[defs->ngenes * (num_nucs - 2) + i] +
               v0[defs->ngenes * (num_nucs - 1) + i] * (-eps / (theta * theta) + 1 / (2.0 * theta) - 2 * delta / theta);
            vdot[defs->ngenes * (num_nucs - 1) + i] /= eps / (theta * theta) + 1 / (2.0 * theta);
         }
      }
   } else {
      error("Dvdt31stJeffNeumann: Bad rule %i sent to Dvdt31stJeffNeumann", rule);
   }
}
