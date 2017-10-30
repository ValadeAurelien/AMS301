#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <stdbool.h>
#include<mpi.h>
#define MIN(a, b) ((a)<(b) ? (a) : (b))


typedef struct ext_grid_t {
    double *grid_read,
	*grid_write,
	*left_write,
        *left_read,
	*top_write,
        *top_read,
	*bottom_write,
        *bottom_read,
	*right_write,
        *right_read,
        *leftR,
        *topR,
        *bottomR,
        *rightR;
} ext_grid_t;

enum border_t { TOP=1, BOTTOM=2, LEFT=4, RIGHT=8, NONE=0 };

unsigned int IND(unsigned nx, unsigned ny, unsigned Ny) {
    return nx*Ny+ny;
}

double isoline(double x, double y)
{
    double theta = 3.141592653589793/2; 
    return sin(theta)*x+cos(theta)*y;
}

double exact_sol(double t, double x, double y,
		double T_max, double kappa, double sigma)
{
    const double x_t = isoline(x, y);
    return T_max/sqrt( 1 + 4*t*kappa/pow(sigma, 2) ) *
	exp( -pow(x_t, 2)/(pow(sigma, 2) + 4*t*kappa) );
}

// remplie notre grille locale à partir de ses coordonnées globales
void exact_sol_grid(ext_grid_t m_grid, unsigned m_Nx, unsigned m_Ny,
		    unsigned i, unsigned j, double dx, double dy,
		    double t, double T_max, double kappa, double sigma)
{
    unsigned nx, ny;
    for (ny=j+1; ny< j+m_Ny-1; ++ny) {
	for (nx=i+1; nx< i+m_Nx-1; ++nx) {
	    m_grid.grid_read[IND(nx-i-1, ny-j-1, m_Ny)] = exact_sol(t, nx*dx, ny*dy,
                                                          T_max, kappa, sigma);
        }
    }
    for(ny=0; ny<m_Ny; ++ny){
        m_grid.bottom_read[ny] = exact_sol(t, (i+m_Nx-1)*dx, (ny+j)*dy, T_max, kappa, sigma);
        m_grid.top_read[ny] = exact_sol(t, i*dx, (ny+j)*dy, T_max, kappa, sigma);
    }
    for(nx=0; nx<m_Nx; ++nx){
        m_grid.left_read[nx] = exact_sol(t, (nx+i)*dx, j*dy, T_max, kappa, sigma);
        m_grid.right_read[nx] = exact_sol(t, (nx+i)*dx, (j+m_Nx-1)*dy, T_max, kappa, sigma);
    }
}

void one_point(ext_grid_t m_grid, unsigned border, 
	       unsigned nx, unsigned ny, unsigned m_Nx, unsigned m_Ny,
	       double dt, double dx, double dy, double kappa)
{
    if((border & TOP) && nx == 0) m_grid.top_write[ny] = m_grid.top_read[ny];
    if((border & BOTTOM) && nx == m_Nx-1) m_grid.bottom_write[ny] = m_grid.bottom_read[ny];
    if((border & LEFT) && ny == 0) m_grid.left_write[nx] = m_grid.left_read[nx];
    if((border & RIGHT) && ny == m_Ny-1) m_grid.right_write[nx] = m_grid.right_read[nx];

    double top=0, bottom=0, left=0, right=0, self=0;
    
    
    if ( nx==0 ) top = m_grid.topR[ny], self = m_grid.top_read[ny];
    else if( nx==1 ) top = m_grid.top_read[ny];
    else if ( nx==m_Nx-1) bottom = m_grid.bottomR[ny], self = m_grid.bottom_read[ny];
    else if( nx==m_Nx-2) bottom = m_grid.bottom_read[ny];
            
    if ( ny==0 ) left = m_grid.leftR[nx], self = m_grid.left_read[nx];
    else if ( ny==1) left = m_grid.left_read[nx];
    else if ( ny==m_Ny-1 ) right = m_grid.rightR[nx], self = m_grid.right_read[nx];
    else if ( ny==m_Ny-2) right = m_grid.right_read[nx];
     
    if(nx > 1) top = m_grid.grid_read[IND(nx-2, ny, m_Ny)];
    if(nx < m_Nx-2) bottom = m_grid.grid_read[IND(nx+1, ny, m_Ny)];
    
    if(ny > 1) left = m_grid.grid_read[IND(nx, ny-2, m_Ny)];
    if(ny < m_Ny-2) right = m_grid.grid_read[IND(nx, ny+1, m_Ny)];
    
    if(nx>0 && nx<m_Nx-1 && ny>0 && ny<m_Ny-1){
        self = m_grid.grid_read[IND(nx, ny, m_Ny)];
        m_grid.grid_write[IND(nx, ny, m_Ny)] = self + dt * kappa *
	( (left+right-2*self)/pow(dx, 2) + (top+bottom-2*self)/pow(dy, 2) ); 
    }else{
        double scheme_calculation = self + dt * kappa *
	( (left+right-2*self)/pow(dx, 2) + (top+bottom-2*self)/pow(dy, 2) );
        if (nx == 0) m_grid.top_write[ny] = scheme_calculation;
        if (nx == m_Nx-1) m_grid.bottom_write[ny] = scheme_calculation;
        if (ny == 0) m_grid.left_write[nx] = scheme_calculation;
        if (ny == m_Ny-1) m_grid.right_write[nx] = scheme_calculation;
    }
    
    return;
}

void one_step(ext_grid_t ext_grid, enum border_t border,
	      unsigned m_Nx, unsigned m_Ny,
	      double dt, double dx, double dy, double kappa)
{
    unsigned nx, ny;
    for (ny=0; ny<m_Ny; ++ny) 
	for (nx=0; nx<m_Nx; ++nx) 
	    one_point(ext_grid, border, 
		      nx, ny, m_Nx, m_Ny,
		      dt, dx, dy, kappa);
}

double L2_diff(const ext_grid_t ext_grid_a, const ext_grid_t ext_grid_b,
	       unsigned m_Nx, unsigned m_Ny)
{
    double diff = 0;
    unsigned nx, ny;
    for (ny=0; ny<m_Ny-2; ++ny) 
	for (nx=0; nx<m_Nx-2; ++nx) 
	    diff += pow(ext_grid_a.grid_read[IND(nx, ny, m_Ny)] -
			ext_grid_b.grid_read[IND(nx, ny, m_Ny)], 2);
    for (ny = 0; ny < m_Ny; ++ny) {
        diff += pow(ext_grid_a.bottom_read[ny] - 
                    ext_grid_b.bottom_read[ny], 2);
        diff += pow(ext_grid_a.top_read[ny] - 
                    ext_grid_b.top_read[ny], 2);
    }
    for (nx = 1; nx < m_Nx-1; ++nx) {
        diff += pow(ext_grid_a.left_read[nx] - 
                    ext_grid_b.left_read[nx], 2);
        diff += pow(ext_grid_a.right_read[nx] - 
                    ext_grid_b.right_read[nx], 2);
    }
    return diff;
}

void init_ext_grid_size(const int threads_nb, int* Nx_threads_nb, 
                        int* Ny_threads_nb, int* m_Nx, int* m_Ny, 
                        int* i, int* j, unsigned int* m_border,
                        int Nx, int Ny, const int rank){
    /*Heuristique sur la dimension des grilles étendues*/
    if((*Nx_threads_nb) * ((*Ny_threads_nb) + 1) <= threads_nb)
        *Ny_threads_nb += 1;
    else if((*Nx_threads_nb + 1) * (*Ny_threads_nb) <= threads_nb)
        *Nx_threads_nb += 1;
    
    *m_Nx = floor(Nx/(*Nx_threads_nb));
    *m_Ny = floor(Ny/(*Ny_threads_nb));
    // Position de notre bout de grille sur la grille réelle
    *i =  (rank/(*Ny_threads_nb)) * (*m_Nx), *j = (rank % (*Ny_threads_nb)) * (*m_Ny);
    //Si l'on se trouve sur le bord droit, on prend toutes les cases de droites restantes
    if((rank + 1) % (*Ny_threads_nb) == 0){
        *m_Ny = Ny - (*m_Ny) * (*Ny_threads_nb - 1);
    }
    // Idem si l'on se retrouve sur le bord du bas, on prend toutes les cases du bas
    if(rank/(*Ny_threads_nb) + 1 == *Nx_threads_nb){
        *m_Nx = Nx - (*m_Nx) * (*Nx_threads_nb - 1);
    }
    *m_border = NONE;
    if(*i == 0) *m_border |= TOP;
    if(*i + *m_Nx == Nx) *m_border |= BOTTOM;
    if(*j == 0) *m_border |= LEFT;
    if(*j + *m_Ny == Ny) *m_border |= RIGHT;
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
    /* Initialisation de MPI */
    
    MPI_Init(&argc, &argv);
    int threads_nb, m_rank, tag = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &threads_nb);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
    
    /* Définition des constantes */
    
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
    double *grid_tmp = NULL, *L2_err = NULL, *times = NULL;
    double L2_err_loc = 0, L2_err_global = 0;
    
    int Nx_threads_nb = floor(sqrt(threads_nb * Nx/Ny)),
        Ny_threads_nb = floor(sqrt(threads_nb * Ny/Nx));
    
    int m_Nx, m_Ny; // Nb de points de discrétisation selon Ox et Oy 
                    // inclus dans notre grille
    int i, j;       // Position de notre grille étendue dans la grille
    unsigned int m_border; 
    init_ext_grid_size(threads_nb, &Nx_threads_nb, &Ny_threads_nb, &m_Nx, &m_Ny, &i, &j, &m_border, Nx, Ny, m_rank);
    
    printf("\n m_border = %u \n m_Nx = %d \n m_Ny = %d \n Nx_threads_nb = %d \n Ny_threads_nb = %d \n i = %d \n j = %d \n",m_border, m_Nx, m_Ny, Nx_threads_nb, Ny_threads_nb, i, j);
    
    /* initialisation de notre bout de grille */
    ext_grid_t m_grid, grid_sol_exact;
    
    // Allocation de mémoire pour le centre et les bords de notre grille
    bool allocating_err = false;
    
    m_grid.grid_read = (double*) malloc((m_Nx - 2) * (m_Ny - 2) * sizeof(double));
    m_grid.grid_write = (double*) malloc((m_Nx - 2) * (m_Ny - 2) * sizeof (double));
    grid_sol_exact.grid_read = (double*) malloc((m_Nx - 2) * (m_Ny - 2) * sizeof(double));
    allocating_err = !m_grid.grid_read || !m_grid.grid_write || !grid_sol_exact.grid_read;

// TODO: pas besoin d'allouer de mémoire à certaines grilles au bord 
// qui n'enverrons jamais leurs bords.    
    
    m_grid.leftR = (double*) malloc(m_Nx * sizeof (double));
    m_grid.left_read = (double*) malloc(m_Nx * sizeof (double));
    m_grid.left_write = (double*) malloc(m_Nx * sizeof (double));
    grid_sol_exact.left_read = (double*) malloc(m_Nx * sizeof (double));
    allocating_err |= !m_grid.leftR || !m_grid.left_read 
                      || !m_grid.left_write || !grid_sol_exact.left_read;

    m_grid.rightR = (double*) malloc(m_Nx * sizeof (double));
    m_grid.right_read = (double*) malloc(m_Nx * sizeof (double));
    m_grid.right_write = (double*) malloc(m_Nx * sizeof (double));
    grid_sol_exact.right_read = (double*) malloc(m_Nx * sizeof (double));
    allocating_err |= !m_grid.rightR || !m_grid.right_read 
                      || !m_grid.right_write || !grid_sol_exact.right_read;

    m_grid.topR = (double*) malloc(m_Ny * sizeof (double));
    m_grid.top_read = (double*) malloc(m_Ny * sizeof (double));
    m_grid.top_write = (double*) malloc(m_Ny * sizeof (double));
    grid_sol_exact.top_read = (double*) malloc(m_Ny * sizeof (double));
    allocating_err |= !m_grid.topR || !m_grid.top_read 
                      || !m_grid.top_write || !grid_sol_exact.top_read;

    m_grid.bottomR = (double*) malloc(m_Ny * sizeof (double));
    m_grid.bottom_read = (double*) malloc(m_Ny * sizeof (double));
    m_grid.bottom_write = (double*) malloc(m_Ny * sizeof (double));
    grid_sol_exact.bottom_read = (double*) malloc(m_Ny * sizeof (double));
    allocating_err |= !m_grid.bottomR || !m_grid.bottom_read 
                      || !m_grid.bottom_write || !grid_sol_exact.bottom_read;
    
    
    // Le thread 0 reçoit les erreurs des autres threads et calcule l'erreur.
    if (calc_err && m_rank == 0) {
	L2_err = (double*) malloc(nb_steps*sizeof(double));
	times = (double*) malloc(nb_steps*sizeof(double));
        allocating_err |= !times || !L2_err;
    } 
    
    if(allocating_err){
        printf("Cannot allocate memory\n");
        return EXIT_FAILURE;
    }
    
    const double dt = t_end/nb_steps,
	dx = a/Nx,
	dy = b/Ny,
	CFL = kappa*dt/MIN(pow(dx, 2), pow(dy, 2));
    if(m_rank == 0) printf("CFL : %f \n", CFL);
    if (CFL>CFL_max) {
	if(m_rank == 0) printf("CFL too big... Aborting! \n");
	return EXIT_FAILURE;
    }
    
    /*
      initialisation des grilles (chaque thread initialise son bout de grille)
    */
    exact_sol_grid(m_grid, m_Nx, m_Ny, i, j, dx, dy,
		   0, T_max, kappa, sigma);
    /*
      boucle principale
    */
    unsigned step;
    for (step=0; step<nb_steps; step++) {

        // Envoi à droite
        if(m_border & RIGHT == 0)
            MPI_Send(m_grid.right_read, m_Nx, MPI_DOUBLE, m_rank+1, tag, MPI_COMM_WORLD);
        if(m_border & LEFT == 0)
            MPI_Recv(m_grid.leftR, m_Nx, MPI_DOUBLE, m_rank-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //Envoi à gauche
        if(m_border & LEFT == 0)
            MPI_Send(m_grid.left_read, m_Nx, MPI_DOUBLE, m_rank-1, tag, MPI_COMM_WORLD);
        if(m_border & RIGHT == 0)
            MPI_Recv(m_grid.rightR, m_Nx, MPI_DOUBLE, m_rank+1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //Envoi en haut
        if(m_border & TOP == 0)
            MPI_Send(m_grid.top_read, m_Ny, MPI_DOUBLE, m_rank - Ny_threads_nb, tag, MPI_COMM_WORLD);
        if(m_border & BOTTOM == 0)
            MPI_Recv(m_grid.bottomR, m_Ny, MPI_DOUBLE, m_rank + Ny_threads_nb, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //Envoi en bas
        if(m_border & BOTTOM == 0)
            MPI_Send(m_grid.bottom_read, m_Ny, MPI_DOUBLE, m_rank + Ny_threads_nb, tag, MPI_COMM_WORLD);
        if(m_border & TOP == 0)
            MPI_Recv(m_grid.topR, m_Ny, MPI_DOUBLE, m_rank - Ny_threads_nb, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        one_step(m_grid, m_border, m_Nx, m_Ny, dt, dx, dy, kappa);   
         // On avance d'un pas dans le temps
	grid_tmp = m_grid.grid_write;
	m_grid.grid_write = m_grid.grid_read;
	m_grid.grid_read = grid_tmp;
        
        grid_tmp = m_grid.right_write;
        m_grid.right_write = m_grid.right_read;
        m_grid.right_read = grid_tmp;
        
        grid_tmp = m_grid.top_write;
        m_grid.top_write = m_grid.top_read;
        m_grid.top_read = grid_tmp;
        
        grid_tmp = m_grid.bottom_write;
        m_grid.bottom_write = m_grid.bottom_read;
        m_grid.bottom_read = grid_tmp;
        
        grid_tmp = m_grid.left_write;
        m_grid.left_write = m_grid.left_read;
        m_grid.left_read = grid_tmp;
        
        // Réception de toutes les erreurs sur le thread 0.
        if(calc_err){
             exact_sol_grid(grid_sol_exact, m_Nx, m_Ny, i, j, dx, dy,
			   step*dt, T_max, kappa, sigma);
             if(m_rank > 0) {
                 // L2_err_loc est simplement la somme des carrés des diff
                 // entre la solution réelle et son approximation sur notre 
                 // grille.
                 L2_err_loc = (m_rank < Nx_threads_nb * Ny_threads_nb ? L2_diff(m_grid, grid_sol_exact, m_Nx, m_Ny) : 0);
                 MPI_Reduce(&L2_err_loc, &L2_err_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
             }
             else {
                 L2_err_loc = L2_diff(m_grid, grid_sol_exact, m_Nx, m_Ny);
                 MPI_Reduce(&L2_err_loc, &L2_err_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                 L2_err[step] = sqrt(L2_err_global)/(Nx*Ny); 
                 times[step] = step*dt;
                 printf("Time : %f \n Error : %f \n", times[step], L2_err[step]);
             }
        }
    }
 /*
    free(grid_sol_exact.bottom_read);
    free(m_grid.bottom_write);
    free(m_grid.bottom_read);
    free(m_grid.bottomR);
    free(grid_sol_exact.top_read);
    free(m_grid.top_write);
    free(m_grid.top_read);
    free(m_grid.topR);
    free(grid_sol_exact.right_read);
    free(m_grid.right_write);
    free(m_grid.right_read);
    free(m_grid.rightR);
    free(grid_sol_exact.left_read);
    free(m_grid.left_write);
    free(m_grid.left_read);
    free(m_grid.leftR);
    free(grid_sol_exact.grid_read);
    free(m_grid.grid_write);
    free(m_grid.grid_read);
    free(times);
    free(L2_err);
*/
    MPI_Finalize();

    return EXIT_SUCCESS;

}
    