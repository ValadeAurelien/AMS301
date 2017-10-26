#include <stdio.h>
#include <stdlib.h>

int plot_curve_1D(FILE* pipe, const double* X, const double* Y,
		  unsigned size, char* pre_commands, char* options)
{
    fprintf(pipe, pre_commands);
    fprintf(pipe, "plot '-' ");
    fprintf(pipe, options);
    fprintf(pipe, ";\n");
    unsigned n;
    for (n=0; n<size; n++) {
	fprintf(pipe, "%g %g\n", (float) X[n], (float) Y[n]);
    }
    fprintf(pipe, "e\n");
    fflush(pipe);
    return EXIT_SUCCESS;
}

/* FILE* gnuplot_pipe = popen("gnuplot -p", "w"); */
int plot_grid_2D_map(FILE* pipe, const double* grid,
	    const unsigned Nx, const unsigned Ny,
	    const double dx, const double dy,
	    const char* pre_commands, const char* options)
{
    fprintf(pipe, pre_commands);
    fprintf(pipe, "set pm3d map; splot '-' ");
    fprintf(pipe, options);
    fprintf(pipe, ";\n");
    unsigned nx, ny;
    for (nx=0; nx<Ny; nx++) {
	for (ny=0; ny<Nx; ny++) {
	    /* printf("%g %g %g\n", */
	    /* 	   (float) nx*dx, (float) ny*dy, (float) grid[ny*Nx+nx]); */
	    fprintf(pipe, "%g %g %g\n",
		    (float) nx*dx, (float) ny*dy, (float) grid[ny*Nx+nx]);
	}
	fprintf(pipe, "\n");
    }
    fprintf(pipe, "e\n");
    fflush(pipe);
    return EXIT_SUCCESS;
}
		    
void show_grid(const double* grid, const unsigned Nx, const unsigned Ny)
{
    unsigned nx, ny;
    for (ny=0; ny<Ny; ny++) {
	for (nx=0; nx<Nx; nx++) {
	    printf("%2.2g   ", grid[ny*Nx+nx]);
	}
	printf("\n");
    }
}
