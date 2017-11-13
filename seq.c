#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "utils.h"
#define MIN(a, b) ((a)<(b) ? (a) : (b))


unsigned to_ind(unsigned nx, unsigned ny, unsigned Nx)
{
    return Nx*ny+nx;
}

double ray2(double x, double y, double a, double b)
{
    return pow(x-a/2., 2)+pow(y-b/2., 2);
}

double isoline(double x, double y)
{
    double theta = 3.141592653589793/2; 
    return sin(theta)*x+cos(theta)*y;
}

double exact_sol(double t, double x, double y,
		double T_max, double kappa, double sigma,
		double a, double b)
{
    const double x_t = isoline(x, y); //ray2(x, y, a, b);
    return T_max/sqrt( 1 + 4*t*kappa/pow(sigma, 2) ) *
	exp( -pow(x_t, 2)/(pow(sigma, 2) + 4*t*kappa) );
}

void exact_sol_grid(double* grid, unsigned Nx, unsigned Ny,
		    double dx, double dy,
		    double t, double T_max, double kappa, double sigma,
		    double a, double b)
{
    unsigned nx, ny;
    for (ny=0; ny<Ny; ++ny) {
	for (nx=0; nx<Nx; ++nx) {
	    grid[to_ind(nx, ny, Nx)] = exact_sol(t, nx*dx, ny*dy,
						 T_max, kappa, sigma,
						 a, b);
	}
    }
}

void one_point(double* grid_write, const double* grid_read,
	       unsigned nx, unsigned ny, unsigned Nx, unsigned Ny,
	       double dt, double dx, double dy, double kappa)
{
    if ( nx==0 || ny==0 || nx==Nx-1 || ny==Ny-1) {
	grid_write[to_ind(nx, ny, Nx)] = grid_read[to_ind(nx, ny, Nx)];
	return;
    }
    else {
	const double left = grid_read[to_ind(nx-1, ny, Nx)],
	    top = grid_read[to_ind(nx, ny-1, Nx)],
	    right = grid_read[to_ind(nx+1, ny, Nx)],
	    bottom = grid_read[to_ind(nx, ny+1, Nx)],
	    self = grid_read[to_ind(nx, ny, Nx)];
	grid_write[to_ind(nx, ny, Nx)] = self + dt * kappa *
	    ( (left+right-2*self)/pow(dx, 2) + (top+bottom-2*self)/pow(dy, 2) );
	return;
    }
}

void one_step(double* grid_write, const double* grid_read,
	      unsigned Nx, unsigned Ny,
	      double dt, double dx, double dy, double kappa)
{
    unsigned nx, ny;
    for (ny=0; ny<Ny; ++ny) {
	for (nx=0; nx<Nx; ++nx) {
	    one_point(grid_write, grid_read,
		      nx, ny, Nx, Ny,
		      dt, dx, dy, kappa);
	}
    }
}

double L2_diff(const double* grid_a, const double* grid_b,
	       unsigned Nx, unsigned Ny)
{
    double diff;
    unsigned nx, ny;
    for (ny=0; ny<Ny; ++ny) {
	for (nx=0; nx<Nx; ++nx) {
	    diff += pow(grid_a[to_ind(nx, ny, Nx)]-grid_b[to_ind(nx, ny, Nx)], 2);
	}
    }
    return sqrt(diff)/(Nx*Ny);
}

int main(int argc, char **argv)
{
    /*
      lecture des arguments
    */
    if (argc<13) {
	printf("Missing args (got %d instead of 12) : \n", argc);
        printf("Args : Nx, Ny, nb_steps, t_end, a, b, kappa, CFL, T_max, sigma, plot_every, calc_err\n");
	return EXIT_FAILURE;
    }
    
    const int Nx = floor(atof(argv[1])),
	Ny = floor(atof(argv[2])),
	nb_steps = floor(atof(argv[3]));
    const double t_end = atof(argv[4]),
	a = atof(argv[5]),
	b = atof(argv[6]),
	kappa = atof(argv[7]),
	CFL_max = atof(argv[8]),
	T_max = atof(argv[9]),
	sigma = atof(argv[10]);
    const int plot_every = atoi(argv[11]),
	calc_err = atoi(argv[12]);
    /*
      initialisation des variables
    */
    double *grid_a, *grid_b,
	*grid_sol_exact,
	*grid_read, *grid_write, *grid_tmp,
	*L2_err, *times;
    grid_a = (double*) malloc(Nx*Ny*sizeof(double));
    grid_b = (double*) malloc(Nx*Ny*sizeof(double));
    grid_sol_exact = (double*) malloc(Nx*Ny*sizeof(double));
    grid_read = grid_a;
    grid_write = grid_b;
    if (!grid_a || !grid_b) {
	printf("Cannot allocate memory\n");
	return EXIT_FAILURE;
    }
    if (calc_err) {
	L2_err = (double*) malloc(nb_steps*sizeof(double));
	times = (double*) malloc(nb_steps*sizeof(double));
    }
    const double dt = t_end/nb_steps,
	dx = a/Nx,
	dy = b/Ny,
	CFL = kappa*dt/MIN(pow(dx, 2), pow(dy, 2));
    printf("CFL : %f \n", CFL);
    if (CFL>CFL_max) {
	printf("CFL too big... Aborting! \n");
	return EXIT_FAILURE;
    }
    FILE* gnuplot_pipe_2D,
	* gnuplot_pipe_err;
    char gpt2D_pre_cmd[256],
	gpterr_pre_cmd[256],
	gpt2D_opt[64] = "t '0'",
	gpterr_opt[64] = "w l lc 6";
    if (plot_every) { 
	gnuplot_pipe_2D = popen("gnuplot -p", "w");
	sprintf(gpt2D_pre_cmd, "set yrange [0:%f];", T_max);
    }    
    if (calc_err) { 
	gnuplot_pipe_err = popen("gnuplot -p", "w");
    }
    
    /*
      initialisation des grilles
    */
    unsigned nx, ny;
    exact_sol_grid(grid_sol_exact, Nx, Ny, dx, dy,
		   0, T_max, kappa, sigma,
		   a, b);
    memcpy(grid_a, grid_sol_exact, Nx*Ny*sizeof(double));
    if (plot_every) {
	plot_grid_2D_map(gnuplot_pipe_2D, grid_read,
			 Nx, Ny, dx, dy, gpt2D_pre_cmd, gpt2D_opt);
    }
    if (calc_err) {
	times[0] = 0;
	L2_err[0] = 0;
    }
    /*
      boucle principale
    */
    unsigned step;
    for (step=0; step<nb_steps; step++) {

	one_step(grid_write, grid_read,
		 Nx, Ny, dt, dx, dy,
		 kappa);
	grid_tmp = grid_write;
	grid_write = grid_read;
	grid_read = grid_tmp;

	if (plot_every) {
	    if (!(step%plot_every)) {
	    	sprintf(gpt2D_opt, "t 't=%f'", step*dt);
		plot_grid_2D_map(gnuplot_pipe_2D, grid_read,
				 Nx, Ny, dx, dy, gpt2D_pre_cmd, gpt2D_opt);
	    }
	}
	if (calc_err) {
	    exact_sol_grid(grid_sol_exact, Nx, Ny, dx, dy,
			   step*dt, T_max, kappa, sigma,
			   a, b);
	    times[step] = step*dt;
	    L2_err[step] = L2_diff(grid_read, grid_sol_exact, Nx, Ny);
	    if ( !(step%(int) floor(nb_steps/1000.)) ) {
	        printf("t = %f , erreur = %f \n", times[step], L2_err[step]); fflush(stdout);
	    }
	}
    }
    double exec_time = (double)clock()/CLOCKS_PER_SEC;
    printf("\n Temps d'exécution du programme : %f s\n", exec_time); fflush(stdout);
    if (calc_err) {
	plot_curve_1D(gnuplot_pipe_err, times, L2_err,
		      nb_steps, gpterr_pre_cmd, gpterr_opt);
    }
    free(grid_a);
    free(grid_b);
    free(grid_sol_exact);
    free(times);
    free(L2_err);
    pclose(gnuplot_pipe_2D);
    pclose(gnuplot_pipe_err);
    return EXIT_SUCCESS;
}
    
    
