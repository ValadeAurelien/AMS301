#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include "utils.h"
#define MIN(a, b) ((a)<(b) ? (a) : (b))
#define IND(nx, ny, Nx) Nx*ny+nx

typedef struct ext_grid_t {
    double *grid_read,
	*grid_write,
	*left,
	*top,
	*bottom,
	*right;
} ext_grid_t;

enum border_t { TOP=1, BOTTOM=2, LEFT=3, RIGHT=4, NONE=0 };

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
    for (ny=0; ny<Ny; ++ny) 
	for (nx=0; nx<Nx; ++nx) 
	    grid[IND(nx, ny, Nx)] = exact_sol(t, nx*dx, ny*dy,
						 T_max, kappa, sigma,
						 a, b);
}

void one_point(ext_grid_t ext_grid, border_t border, 
	       unsigned nx, unsigned ny, unsigned Nx, unsigned Ny,
	       double dt, double dx, double dy, double kappa)
{
    switch ( border ) {
    case TOP:
	if ( ny==0 ) {
	    grid_write[IND(nx, ny, Nx)] = grid_read[IND(nx, ny, Nx)];
	    return;
	}
	break;
    case BOTTOM:
	if ( ny==Ny-1 ) {
	    grid_write[IND(nx, ny, Nx)] = grid_read[IND(nx, ny, Nx)];
	return;
	}
	break;
    case LEFT:
	if ( nx==0 ) {
	    grid_write[IND(nx, ny, Nx)] = grid_read[IND(nx, ny, Nx)];
	return;
	}
	break;
    case RIGHT:
	if ( nx==Nx-1 ) {
	    grid_write[IND(nx, ny, Nx)] = grid_read[IND(nx, ny, Nx)];
	return;
	}
	break;
    case NONE;
    default:
	break;
    }
    const double top,
	bottom,
	left,
	right,
	self;
    if ( ny==0 ) top = ext_grid.top[nx];
    else top = ext_grid.grid_read[IND(nx, ny-1, Nx)];

    if ( ny==Ny-1 ) bottom = ext_grid.bottom[nx];
    else bottom = ext_grid.grid_read[IND(nx, ny+1, Nx)];

    if ( nx==0 ) left = ext_grid.left[ny];
    else left = ext_grid.grid_read[IND(nx-1, ny, Nx)];

    if ( nx==Nx-1 ) right = ext_grid.right[ny];
    else right = ext_grid.grid_read[IND(nx+1, ny, Nx)];

    self = ext_grid.grid_read[IND(nx, ny, Nx)];
    grid_write[IND(nx, ny, Nx)] = self + dt * kappa *
	( (left+right-2*self)/pow(dx, 2) + (top+bottom-2*self)/pow(dy, 2) );
    return;
}

void one_step(ext_grid_t ext_grid, border_t border,
	      unsigned Nx, unsigned Ny,
	      double dt, double dx, double dy, double kappa)
{
    unsigned nx, ny;
    for (ny=0; ny<Ny; ++ny) 
	for (nx=0; nx<Nx; ++nx) 
	    one_point(ext_grid, border, 
		      nx, ny, Nx, Ny,
		      dt, dx, dy, kappa);
}

double L2_diff(const ext_grid_t ext_grid_a, const ext_grid_t ext_grid_b,
	       unsigned Nx, unsigned Ny)
{
    double diff;
    unsigned nx, ny;
    for (ny=0; ny<Ny; ++ny) 
	for (nx=0; nx<Nx; ++nx) 
	    diff += pow(ext_grid_a.grid_read[IND(nx, ny, Nx)] -
			ext_grid_b.grid_read[IND(nx, ny, Nx)], 2);
    return sqrt(diff)/(Nx*Ny);
}

int main(int argc, char **argv)
{
    /*
      lecture des arguments
    */
    if (argc<13) {
	printf("Missing args (got %d instead of 12) : \n");
        printf("Args : Nx, Ny, nb_steps, t_end, a, b, kappa, CFL, T_max, sigma, plot_every, calc_err\n", argc);
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
    exact_sol_grid(m_grid.grid_read, Nx, Ny, dx, dy,
		   0, T_max, kappa, sigma,
		   a, b);
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

	if (calc_err) {
	    exact_sol_grid(grid_sol_exact, Nx, Ny, dx, dy,
			   step*dt, T_max, kappa, sigma,
			   a, b);
	    times[step] = step*dt;
	    L2_err[step] = L2_diff(grid_read, grid_sol_exact, Nx, Ny);
	}
    }

    free(grid_a);
    free(grid_b);
    free(m_grid.left);
    free(m_grid.right);
    free(m_grid.top);
    free(m_grid.bottom);
    free(times);
    free(L2_err);
    return EXIT_SUCCESS;
}
    
    
