#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif

#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif

#define LOWER_TAIL_ONE 7.5
#define UPPER_TAIL_ZERO 20
#ifndef ZEPS
#define ZEPS 1e-20
#endif

#ifndef square
#define square(x) ((x) * (x))
#endif

double ninv(double p)
/****************************************************
   C Equivalent of Wichura's PPND16, Algorithm AS241
   Applied Statistics Vol 37 1988 pp 477 - 484
*****************************************************/
{
   const double SPLIT1 = 0.425,
                SPLIT2 = 5.0,
                CONST1 = 0.180625,
                CONST2 = 1.6;

   static const double a[8] = {
       3.3871328727963666080E0,
       1.3314166789178437745E2,
       1.9715909503065514427E3,
       1.3731693765509461125E4,
       4.5921953931549871457E4,
       6.7265770927008700853E4,
       3.3430575583588128105E4,
       2.5090809287301226727E3};

   static const double b[7] = {
       4.2313330701600911252E1,
       6.8718700749205790830E2,
       5.3941960214247511077E3,
       2.1213794301586595867E4,
       3.9307895800092710610E4,
       2.8729085735721942674E4,
       5.2264952788528545610E3};

   static const double c[8] = {
       1.42343711074968357734E0,
       4.63033784615654529590E0,
       5.76949722146069140550E0,
       3.64784832476320460504E0,
       1.27045825245236838258E0,
       2.41780725177450611770E-1,
       2.27238449892691845833E-2,
       7.74545014278341407640E-4};

   static const double d[7] = {
       2.05319162663775882187E0,
       1.67638483018380384940E0,
       6.89767334985100004550E-1,
       1.48103976427480074590E-1,
       1.51986665636164571966E-2,
       5.47593808499534494600E-4,
       1.05075007164441684324E-9};

   static const double e[8] = {
       6.65790464350110377720E0,
       5.46378491116411436990E0,
       1.78482653991729133580E0,
       2.96560571828504891230E-1,
       2.65321895265761230930E-2,
       1.24266094738807843860E-3,
       2.71155556874348757815E-5,
       2.01033439929228813265E-7};

   static const double f[7] = {
       5.99832206555887937690E-1,
       1.36929880922735805310E-1,
       1.48753612908506148525E-2,
       7.86869131145613259100E-4,
       1.84631831751005468180E-5,
       1.42151175831644588870E-7,
       2.04426310338993978564E-15};

   double q = p - 0.5;
   double r, x;

   if (fabs(q) < SPLIT1)
   {
      r = CONST1 - q * q;
      return q * (((((((a[7] * r + a[6]) * r + a[5]) * r + a[4]) * r + a[3]) * r + a[2]) * r + a[1]) * r + a[0]) /
             (((((((b[6] * r + b[5]) * r + b[4]) * r + b[3]) * r + b[2]) * r + b[1]) * r + b[0]) * r + 1.0);
   }
   else
   {
      if (q < 0)
         r = p;
      else
         r = 1.0 - p;

      if (r < 1e-10)
         return (q < 0 ? -20.0 : 20.0);

      if (r > 0.0)
      {
         r = sqrt(-log(r));
         if (r <= SPLIT2)
         {
            r -= CONST2;
            x = (((((((c[7] * r + c[6]) * r + c[5]) * r + c[4]) * r + c[3]) * r + c[2]) * r + c[1]) * r + c[0]) /
                (((((((d[6] * r + d[5]) * r + d[4]) * r + d[3]) * r + d[2]) * r + d[1]) * r + d[0]) * r + 1.0);
         }
         else
         {
            r -= SPLIT2;
            x = (((((((e[7] * r + e[6]) * r + e[5]) * r + e[4]) * r + e[3]) * r + e[2]) * r + e[1]) * r + e[0]) /
                (((((((f[6] * r + f[5]) * r + f[4]) * r + f[3]) * r + f[2]) * r + f[1]) * r + f[0]) * r + 1.0);
         }
      }
      else
         x = 9;

      if (q < 0)
         x = -x;
      return x;
   }
}
/////////////////////////////////////////
double update_integral(double (*f)(double x, void *extras), void *extras,
                       double a, double b, double previous, int round)
{
   double h, sum;
   int n = 1 << (round - 1);

   if (round == 0)
      return 0.5 * ((*f)(a, extras) + (*f)(b, extras)) * (b - a);

   sum = previous * n / (b - a);
   h = (b - a) / (2 * n);
   int i = 1;
   for (i = 1; i < 2 * n; i += 2)
   {
      sum += (*f)(a + i * h, extras);
   }
   return sum / (2 * n) * (b - a);
}
///////////////////////////////
double integral(double (*f)(double x, void *extras), void *extras,
                double a, double b, double eps, int minrounds)
{
   double old = update_integral(f, extras, a, b, 0.0, 0), result;
   int round = 1;

   while (1)
   {
      result = update_integral(f, extras, a, b, old, round++);
      if (round > minrounds && fabs(result - old) < eps * (fabs(result) + fabs(old)) + ZEPS)
         return result;
      old = result;
   }
}

/////////////////////////////////
double ndist(double z, int upper)
{
   // C version of ID Hill, "The Normal Integral"
   // Applied Statistics, Vol 22, pp. 424-427

   // If 7 digit accuracy is enough, alternative is
   // return erfcc(x / M_SQRT2) * 0.5;

   if (z < 0)
   {
      upper = !upper;
      z = -z;
   }

   if ((z > LOWER_TAIL_ONE && !upper) || z > UPPER_TAIL_ZERO)
      return (upper) ? 0.0 : 1.0;

   double p, y = 0.5 * z * z;

   if (z < 1.28)
   {
      p = 0.5 - z * (0.398942280444 - 0.399903438504 * y /
                                          (y + 5.75885480458 - 29.8213557808 / (y + 2.62433121679 + 48.6959930692 / (y + 5.92885724438))));
   }
   else
   {
      p = 0.398942270385 * exp(-y) /
          (z - 2.8052e-8 + 1.00000615302 / (z + 3.98064794e-4 + 1.98615381364 / (z - 0.151679116635 + 5.29330324926 / (z + 4.8385912808 - 15.1508972451 / (z + 0.742380924027 + 30.789933034 / (z + 3.99019417011))))));
   }

   return (upper) ? p : 1 - p;
}

/////////////////////////////////////////
double stuff_to_integrate(double z1, void *parameters)
{
   double *vector = (double *)parameters;

   double pisamples = vector[0];
   double ncp1 = vector[1];
   double ncp2 = vector[2];
   // double c1 = vector[3];
   double c2 = vector[4];

   double mu = ncp2 * sqrt(1.0 - pisamples) + z1 * sqrt(pisamples);
   double var = 1.0 - pisamples;

   double density = 1.0 / sqrt(2 * M_PI) * exp(-0.5 * (z1 - ncp1) * (z1 - ncp1));

   double z_high = (c2 - mu) / sqrt(var);
   double z_low = (-c2 - mu) / sqrt(var);

   return density * (ndist(z_high, 1) + ndist(z_low, 0));
}

// The main parameter is:
//    cjoint -- the threshold for joint analysis
//
// Additional parameters in the parameter vector are:
//   pisamples = parameters[0]    -- the proportion of samples in stage 1
//   ncp1 = parameters[1]         -- the expected value of the stage 1 statistic
//   ncp2 = parameters[2]         -- the expected value of the stage 2 statistic
//   c1 = parameters[3]           -- threshold for stage 1 analysis
//   c2 = parameters[4]           -- threshold for stage 2 analysis
//

/////////////////////////////////////////
double pjoint(double cjoint, void *parameters)
{
   double *vector = (double *)parameters;

   // double pisamples = vector[0];
   // double ncp1 = vector[1];
   // double ncp2 = vector[2];
   double c1 = vector[3];

   vector[4] = cjoint;

   return (integral(stuff_to_integrate, (void *)vector, c1, c1 + max(c1 + 5, cjoint + 3), 1e-7, 5) +
           integral(stuff_to_integrate, (void *)vector, min(-c1 - 5, -cjoint - 3), -c1, 1e-7, 5));
}

double solve(double (*func)(double, void *extras), void *extras,
             double target, double lo, double hi, double e)
{
   if ((*func)(lo, extras) > (*func)(hi, extras))
   {
      double temp = lo;
      lo = hi;
      hi = temp;
   }

   e = e * target + ZEPS;

   while (1)
   {
      double d = hi - lo;
      double point = lo + d * 0.5;
      double fpoint = (*func)(point, extras);

      if (fpoint < target)
      {
         d = lo - point;
         lo = point;
      }
      else
      {
         d = point - hi;
         hi = point;
      }

      if (fabs(d) < e || fpoint == target)
         return point;
   }
}
