#include <R.h>
#include <Rinternals.h>
//#include <R_ext/Rdynload.h>
#include "cats.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

SEXP cats_(SEXP Rfreq, SEXP Rfreq2, SEXP Rncases, SEXP Rncontrols, SEXP Rncases2, SEXP Rncontrols2, SEXP Rrisk, SEXP Rrisk2, SEXP Rpisamples, SEXP Rprevalence, SEXP Rprevalence2, SEXP Radditive, SEXP Rrecessive, SEXP Rdominant, SEXP Rmultiplicative, SEXP Ralpha, SEXP Rpimarkers)
{

  double freq = REAL(Rfreq)[0];

  int ncases = INTEGER(Rncases)[0];
  int ncontrols = INTEGER(Rncontrols)[0];
  int ncases2 = INTEGER(Rncases2)[0];
  int ncontrols2 = INTEGER(Rncontrols2)[0];

  double pisamples = REAL(Rpisamples)[0];
  double pimarkers = REAL(Rpimarkers)[0];

  double freq2 = REAL(Rfreq2)[0];
  double risk = REAL(Rrisk)[0];
  double risk2 = REAL(Rrisk2)[0];
  double prevalence = REAL(Rprevalence)[0];
  double prevalence2 = REAL(Rprevalence2)[0];

  double alpha = REAL(Ralpha)[0];
  int additive = INTEGER(Radditive)[0];
  int dominant = INTEGER(Rdominant)[0];
  int recessive = INTEGER(Rrecessive)[0];
  int multiplicative = INTEGER(Rmultiplicative)[0];

  if (freq2 == -1)
    freq2 = freq;
  if (risk2 == -1)
    risk2 = risk;
  if (prevalence2 == -1)
    prevalence2 = prevalence;
  double C = -ninv(alpha * 0.5);
  double C1 = -ninv(pimarkers * 0.5);
  double Crep = -ninv(alpha / pimarkers);

  // Setup parameter array for our fancy functions ...
  double parameters[6];

  parameters[0] = 1.0 * (ncases + ncontrols) / (ncases + ncontrols + ncases2 + ncontrols2);
  parameters[1] = parameters[2] = 0.0; // Expected stage 1 and 2 statistics
  parameters[3] = C1;
  parameters[4] = 0.0;
  double Cjoint = solve(pjoint, parameters, alpha, Crep, C, 1e-6);

  double p[3];  // The three genotype frequencies
  double p2[3]; // The three genotype frequencies
  double f[3];  // The three penetrances
  double f2[3]; // The three penetrances

  // Genotype frequencies
  p[0] = square(freq);
  p[1] = 2 * freq * (1. - freq);
  p[2] = square(1. - freq);
  p2[0] = square(freq2);
  p2[1] = 2 * freq2 * (1. - freq2);
  p2[2] = square(1. - freq2);

  // Relative penetrances for low risk genotype is 1.0
  f[2] = 1.0;
  // Relative penetrances for low risk genotype is 1.0
  f2[2] = 1.0;

  // Check that at least one parameter is enabled
  if (!(additive == 1 || dominant == 1 || multiplicative == 1 || recessive == 1))
  {
    error("No genetic model selected!\n\n"
          "No command line parameter required for the default multiplicative model.\n"
          "Use the --additive, --recessive or --dominant parameters only once to\n"
          "select a different model.\n\n");
  }

  /////////////////stage 1
  // Relative penetrances for high risk genotypes depend on model
  if (additive == 1)
  {
    f[0] = 2.0 * risk - 1.0;
    f[1] = risk;
  }
  if (dominant == 1)
  {
    f[0] = risk;
    f[1] = risk;
  }
  if (recessive == 1)
  {
    f[0] = risk;
    f[1] = 1.0;
  }
  if (multiplicative == 1)
  {
    f[0] = risk * risk;
    f[1] = risk;
  }

  double scale = prevalence / (f[0] * p[0] + f[1] * p[1] + f[2] * p[2]);

  // Adjusted penetrances
  f[0] *= scale;
  f[1] *= scale;
  f[2] *= scale;

  if (f[0] >= 1.0 || f[0] <= 0)
  {
    error("I don't like the genetic model you requested!\n");
  }

  double pcases = (f[0] * p[0] + f[1] * p[1] * 0.5) / prevalence;
  double pcontrols = ((1. - f[0]) * p[0] + (1. - f[1]) * p[1] * 0.5) / (1. - prevalence);

  double vcases = pcases * (1.0 - pcases);
  double vcontrols = pcontrols * (1.0 - pcontrols);
  ////////////stage 2
  // Relative penetrances for high risk genotypes depend on model
  if (additive)
  {
    f2[0] = 2.0 * risk2 - 1.0;
    f2[1] = risk2;
  }
  if (dominant)
  {
    f2[0] = risk2;
    f2[1] = risk2;
  }
  if (recessive)
  {
    f2[0] = risk2;
    f2[1] = 1.0;
  }
  if (multiplicative)
  {
    f2[0] = risk2 * risk2;
    f2[1] = risk2;
  }

  double scale2 = prevalence / (f2[0] * p2[0] + f2[1] * p2[1] + f2[2] * p2[2]);

  // Adjusted penetrances
  f2[0] *= scale2;
  f2[1] *= scale2;
  f2[2] *= scale2;

  if (f2[0] > 1.0)
  {
    error("I don't like the genetic model you requested for stage 2!\n");
  }

  double pcases2 = (f2[0] * p2[0] + f2[1] * p2[1] * 0.5) / prevalence2;
  double pcontrols2 = ((1. - f2[0]) * p2[0] + (1. - f2[1]) * p2[1] * 0.5) / (1. - prevalence2);

  double vcases2 = pcases2 * (1.0 - pcases2);
  double vcontrols2 = pcontrols2 * (1.0 - pcontrols2);
  double ncp1 = (pcases - pcontrols) / sqrt((vcases / ncases + vcontrols / ncontrols) / (2.0));
  double ncp2 = (pcases2 - pcontrols2) / sqrt((vcases2 / ncases2 + vcontrols2 / ncontrols2) / (2.0));
  double ncp = (pcases - pcontrols) / sqrt((vcases / (ncases + ncases2) + vcontrols / (ncontrols + ncontrols2)) * 0.5);

  double P = ndist(-C - ncp, 0) + ndist(C - ncp, 1);
  double P1 = ndist(-C1 - ncp1, 0) + ndist(C1 - ncp1, 1);
  double Prep = ndist(-C1 - ncp1, 0) * ndist(-Crep - ncp2, 0) +
                ndist(C1 - ncp1, 1) * ndist(Crep - ncp2, 1);

  // double var1= ( (vcases  + vcontrols )/( ncases+ ncontrols) / (2.0));
  //double var2= ( (vcases2 + vcontrols2) /( ncases2+ ncontrols2) / ( 2.0));
  double var1 = ((vcases / ncases + vcontrols / ncontrols) / (2.0));
  double var2 = ((vcases2 / ncases2 + vcontrols2 / ncontrols2) / (2.0));

  // Setup the parameters for our weird integral ...
  //parameters[0] = 1.0*(ncases+ncontrols)/(ncases+ncontrols+ncases2+ncontrols2);
  // parameters[0] = (((vcases2+ vcontrols2) /(ncases2+ ncontrols2))/((vcases+vcontrols) / (ncases+ncontrols)+(vcases2+ vcontrols2)/(ncontrols2+ ncontrols2)));
  if (pisamples == -1)
    parameters[0] = 1.0 / var1 / (1.0 / var1 + 1.0 / var2);
  else
    parameters[0] = pisamples;

  parameters[1] = ncp1;
  parameters[2] = ncp2;
  parameters[3] = C1;
  parameters[4] = Cjoint;

  double upper_from = max(C1, ncp1 - 7);
  double upper_to = max(C1, ncp1) + 5;

  double Pjoint = integral(stuff_to_integrate, (void *)parameters, upper_from, upper_to, 1e-7, 5) +
                  integral(stuff_to_integrate, (void *)parameters, -C1 - 5, -C1, 1e-7, 5);
  double min5 = 1 - pow(1 - Pjoint, 5.0);

  SEXP Rpjoint;
  PROTECT(Rpjoint = allocMatrix(REALSXP, 14, 1));
  REAL(Rpjoint)
  [0] = P;
  REAL(Rpjoint)
  [1] = P1;
  REAL(Rpjoint)
  [2] = Prep;
  REAL(Rpjoint)
  [3] = min5;
  REAL(Rpjoint)
  [4] = Pjoint;
  REAL(Rpjoint)
  [5] = parameters[0];
  REAL(Rpjoint)
  [6] = C;
  REAL(Rpjoint)
  [7] = C1;
  REAL(Rpjoint)
  [8] = Crep;
  REAL(Rpjoint)
  [9] = Cjoint;
  REAL(Rpjoint)
  [10] = pcases;
  REAL(Rpjoint)
  [11] = pcontrols;
  REAL(Rpjoint)
  [12] = pcases2;
  REAL(Rpjoint)
  [13] = pcontrols2;
  UNPROTECT(1);

  return (Rpjoint);
}
