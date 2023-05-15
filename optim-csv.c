#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define RND drand48()
//#define JUST_TEST 1   // uncomment to run test with Belen's parameters

#define Inf 99999       // silly large number for infinity

// global parameters
double SCALE;    // width of perturbation kernel
int FORCESYMM;   // whether to enforce symmetry (not used)

// structure for a model parameterisation
typedef struct
{
  double b[10];       // slopes
  double tau[10];     // changepoints
  double c, sigma;    // intercept and noise width
} Params;

// produce gaussian random number
double gsl_ran_gaussian(const double sigma)
{
  double x, y, r2;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      x = -1 + 2 * RND;
      y = -1 + 2 * RND;

      /* see if it is in the unit circle */
      r2 = x * x + y * y;
    }
  while (r2 > 1.0 || r2 == 0);

  /* Box-Muller transform */
  return sigma * y * sqrt (-2.0 * log (r2) / r2);
}

      
// simply read out mean prediction at x from changepoint model
double modelval(Params P, int model, double x)
{
  // general formula is
  // E(y) = c + sum_{i: tau(i) < x} b[i]*(tau[i]-tau[i-1]) + b[j]*(x-tau[j])
  // where j is the phase in which x falls and i is the set of previous phases

  if(x < P.tau[1] || model <= 1) { return P.b[1]*x+P.c; }
  
  else if(x < P.tau[2] || model <= 2) { return P.b[2]*(x-P.tau[1]) + (P.b[1]*P.tau[1] +P.c); }

  else if(x < P.tau[3] || model <= 3) { return P.b[3]*(x-P.tau[2]) + P.b[2]*(P.tau[2]-P.tau[1]) + (P.b[1]*P.tau[1] +P.c); }

  else if(x < P.tau[4] || model <= 4) { return P.b[4]*(x-P.tau[3]) + P.b[3]*(P.tau[3]-P.tau[2]) + P.b[2]*(P.tau[2]-P.tau[1]) + (P.b[1]*P.tau[1] +P.c); }

  else if(x < P.tau[5] || model <= 5) { return P.b[5]*(x-P.tau[4]) + P.b[4]*(P.tau[4]-P.tau[3]) + P.b[3]*(P.tau[3]-P.tau[2]) + P.b[2]*(P.tau[2]-P.tau[1]) + (P.b[1]*P.tau[1] +P.c); }

  else if(x < P.tau[6] || model <= 6) { return P.b[6]*(x-P.tau[5]) + P.b[5]*(P.tau[5]-P.tau[4]) + P.b[4]*(P.tau[4]-P.tau[3]) + P.b[3]*(P.tau[3]-P.tau[2]) + P.b[2]*(P.tau[2]-P.tau[1]) + (P.b[1]*P.tau[1] +P.c); }

  else if(x < P.tau[7] || model <= 7) { return P.b[7]*(x-P.tau[6]) + P.b[6]*(P.tau[6]-P.tau[5]) + P.b[5]*(P.tau[5]-P.tau[4]) + P.b[4]*(P.tau[4]-P.tau[3]) + P.b[3]*(P.tau[3]-P.tau[2]) + P.b[2]*(P.tau[2]-P.tau[1]) + (P.b[1]*P.tau[1] +P.c); }

  else              { return P.b[8]*(x-P.tau[7]) + P.b[7]*(P.tau[7]-P.tau[6]) + P.b[6]*(P.tau[6]-P.tau[5]) + P.b[5]*(P.tau[5]-P.tau[4]) + P.b[4]*(P.tau[4]-P.tau[3]) + P.b[3]*(P.tau[3]-P.tau[2]) + P.b[2]*(P.tau[2]-P.tau[1]) + (P.b[1]*P.tau[1] +P.c); }
}

// log likelihood for a parameterisation given data x, y
double LogLik(Params P, int model, double *x, double *y, int n)
{
  int i;
  double loglik = 0;
  double predicty;

  // lazy sum over the n datapoints
  // in reality as many datapoints have the same x we don't need to compute model predictions each time
  for(i = 0; i < n; i++)
    {
      predicty = modelval(P, model, x[i]);
      loglik += log(1./sqrt(2.*3.14159*P.sigma*P.sigma)) - (y[i]-predicty)*(y[i]-predicty)/(2*P.sigma*P.sigma);
    }
  return loglik;
}

// specific test parameterisations to compare output from a different solver
void BelenP(Params *P, int expt)
{
  if(expt == 0)
    {
      P->sigma = 0.553; P->c = 0.322;
      P->b[1] = 0.246; P->b[2] = -0.023; P->b[3] = 0.0618; P->b[4] = -0.023;
      P->tau[1] = 11.085; P->tau[2] = 92.568; P->tau[3] = 119.512;
    }
  else if(expt == 1)
    {
      P->sigma = 0.626; P->c = 0.8346;
      P->b[1] = 0.145; P->b[2] = -0.034; P->b[3] = 0.01589; P->b[4] = -0.022;
      P->tau[1] = 15.83; P->tau[2] = 65.07; P->tau[3] = 130.1;
    }
  else if(expt == 2)
    {
      P->sigma = 0.723; P->c = 0.43;
      P->b[1] = 0.1528; P->b[2] = -0.0425; P->b[3] = 0.0983; P->b[4] = -0.02;
      P->tau[1] = 13.79; P->tau[2] = 84.982; P->tau[3] = 109.6;
    }
}

// initial parameterisation with which to begin the optimiser
void defaultP(Params *P, int model, double maxx, int CUSTOM_SPACING)
{
  int i;

  // flat lines, evenly spaced changepoints, sigma = 1
  P->sigma = 1;
  for(i = 0; i < 10; i++)
    P->b[i] = 0;
  P->c = 0;
  for(i = 0; i < 10; i++)
    P->tau[i] = maxx*i/model;

  // we can impose manual spacing of changepoints to better match the switches between behaviours which we anticipate
  if(CUSTOM_SPACING)
    {
      if(model == 2) P->tau[1] = 200/2; 
      if(model == 3) { P->tau[1] = 50/2; P->tau[2] = 200/2; }
      if(model == 4) { P->tau[1] = 50/2; P->tau[2] = 200/2; P->tau[3] = 250/2; }
      if(model == 5) { P->tau[1] = 20/2; P->tau[2] = 50/2; P->tau[3] = 200/2; P->tau[4] = 250/2; }
      if(model == 6) { P->tau[1] = 20/2; P->tau[2] = 50/2; P->tau[3] = 200/2; P->tau[4] = 220/2; P->tau[5] = 250/2; }
      if(model == 7) { P->tau[1] = 20/2; P->tau[2] = 50/2; P->tau[3] = 100/2; P->tau[4] = 200/2; P->tau[5] = 220/2; P->tau[6] = 250/2; }
      if(model == 8) { P->tau[1] = 10/2; P->tau[2] = 20/2; P->tau[3] = 50/2; P->tau[4] = 200/2; P->tau[5] = 210/2; P->tau[6] = 250/2; P->tau[7] = 300/2;}
    }
}

// read data from file
// NO ERROR CHECKING: file is assumed to exist and to contain two space-separated columns of numbers
void ReadData(char *fname, double *x, double *y, int *n)
{
  FILE *fp;

  (*n) = 0;
  fp = fopen(fname, "r");
  while(!feof(fp))
    {
      fscanf(fp, "%lf %lf", &(x[*n]), &(y[*n]));
      if(feof(fp)) break;
      (*n)++;
    }
  fclose(fp);
  printf("Read %i\n", *n);
}

// apply a perturbation to a parameter set, to obtain a new trial parameterisation
void PerturbP(Params *P)
{
  int i;

  // the widths of normal kernels applied to each parameter
  double sigmascale = 0.1*SCALE;
  double cscale = 0.05*SCALE;
  double bscale = 0.005*SCALE;
  double tauscale = 5*SCALE;

  // apply those kernels
  P->sigma += gsl_ran_gaussian(sigmascale);
  P->c += gsl_ran_gaussian(cscale);
  if(P->sigma < 0) P->sigma = 1./Inf;   // sigma should be nonzero positive
  for(i = 0; i < 10; i++)
    P->b[i] += gsl_ran_gaussian(bscale);
  for(i = 0; i < 10; i++)
    {
      P->tau[i] += gsl_ran_gaussian(tauscale);
      if(P->tau[i] < P->tau[i-1]) P->tau[i] = P->tau[i-1]; // tau[i] should not come before tau[i-1]
    }
}
  
int main(int argc, char *argv[])
{
  Params P, newP;
  int n;
  double fx[10000], fy[10000];
  double x[10000], y[10000];
  double loglik, newloglik;
  int i;
  FILE *fp, *fp1;
  int MODEL;
  int MODELSET[] = {1, 2, 4, 6, 8};
  int modelindex;
  char fstr[1000];
  int r;
  int boot;
  double temperature;
  double maxx;
  
  srand48(12);

  // process command line arguments
  if(argc != 4)
    {
      printf("Need datafile, scale (e.g. 0.1), force symmetry (0 or 1)\n");
      return 0;
    }

  SCALE = atof(argv[2]);
  FORCESYMM = atoi(argv[3]);

  printf("Working with %s %f %i\n", argv[1], SCALE, FORCESYMM);

  // read data from file
  ReadData(argv[1], fx, fy, &n);

#ifdef JUST_TEST
  // testing block to resolve differences with Belen's modelling
  printf("Testing outputs...\n");
  MODEL = 4;
  BelenP(&P, 0);
  printf("%f\n", LogLik(P, MODEL, fx, fy, n));
  BelenP(&P, 2);
  printf("%f\n", LogLik(P, MODEL, fx, fy, n));
  fp = fopen("belen-test.txt", "w");
  for(i = 0; i < 250; i++)
    fprintf(fp, "%i %f\n", i, modelval(P, 4, i));
  fclose(fp);
  return 0;
#endif
    
  // open summary file for output
  sprintf(fstr, "summary-%s.csv", argv[1]);
  fp1 = fopen(fstr, "w");
  fprintf(fp1, "scale,forcesymm,model,boot.label,log.lik,sigma,c,");
  for(i = 0; i < 10; i++) fprintf(fp1, "b%i,", i);
  for(i = 0; i < 9; i++) fprintf(fp1, "tau%i,", i);
  fprintf(fp1, "tau9\n");
  
  sprintf(fstr, "output-%s.csv", argv[1]);
  fp = fopen(fstr, "w");
  fprintf(fp, "scale,forcesymm,model,boot.label,t,y\n");
  
  // loop through different model structures (number of phases)
  for(modelindex = 0; modelindex < 5; modelindex++)
    {
      MODEL = MODELSET[modelindex];
      
      // open model-specific file for output
     
      // loop over bootstrap resamples
      for(boot = 0; boot < 300; boot++)
	{
	  // construct resampled dataset
	  // for the first instance, retrieve the original data
	  maxx = 0;
	  for(i = 0; i < n; i++)
	    {
	      r = RND*n;
	      x[i] = (boot == 0 ? fx[i] : fx[r]);
	      y[i] = (boot == 0 ? fy[i] : fy[r]);
	      if(x[i] > maxx) maxx = x[i];
	    }

	  // set initial parameters and calculate log likelihood
	  defaultP(&P, MODEL, maxx, 1);
	  loglik = LogLik(P, MODEL, x, y, n);

	  // initialise temperature and counter
	  temperature = 10;
	  i = 0;

	  // annealing loop
	  while(temperature > 1e-7)
	    {
	      // periodic output
	      i++;
	      if(i % 10000 == 0)
		printf("%i %i %i %f %f\n", MODEL, boot, i, temperature, loglik);
	      // try new parameterisation
	      newP = P;
	      PerturbP(&newP);
	      if(FORCESYMM)
		{
		  switch(MODEL)
		    {
		    case 4: P.b[3] = P.b[1]; P.b[4] = P.b[2]; break;
		    case 6: P.b[4] = P.b[1]; P.b[5] = P.b[2]; P.b[6] = P.b[3]; break;
		    case 8: P.b[5] = P.b[1]; P.b[6] = P.b[2]; P.b[7] = P.b[3]; P.b[8] = P.b[4]; break;
		    }
		}
	      newloglik = LogLik(newP, MODEL, x, y, n);
	      // metropolis acceptance of new parameterisation
	      if(newloglik > loglik || exp(-(loglik-newloglik)/temperature) > RND)
		{
		  loglik = newloglik;
		  P = newP;
		}
	      // cool temperature
	      temperature /= 1.0001;
	      //      printf("%i %f : %f %f\n", MODEL, loglik, P.sigma, P.b[1]);
	    }

	  // output final parameterisation for this resample
	  fprintf(fp1, "%.3f,%i,%i,%i,%f,%f,%f,", SCALE, FORCESYMM, MODEL, boot, loglik, P.sigma, P.c);
	  for(i = 0; i < 10; i++)
	    fprintf(fp1, "%f,", P.b[i]);
	  for(i = 0; i < 9; i++)
	    fprintf(fp1, "%f,", P.tau[i]);
	  fprintf(fp1, "%f\n", P.tau[9]);

	  // output explicit time series for final parameterisation for this bootstrap
	  for(i = 0; i < 250; i++)
	    fprintf(fp, "%.3f,%i,%i,%i,%i,%f\n", SCALE, FORCESYMM, MODEL, boot, i, modelval(P, MODEL, i));
	  fprintf(fp, "\n");
	}
    }
  fclose(fp1);
  fclose(fp);
  
  return 0;
}
  
  
