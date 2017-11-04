#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <stdbool.h>
#include <time.h>
#include <mpi.h>
#include "utils.h"
#define MIN(a, b) ((a)<(b) ? (a) : (b))
/* 
   CODES ERREURS 
   0 -> no error
   1 -> missing args 
   2 -> cannot allocate mem
   3 -> CFL too big
*/
/* Définition et fonctions sur les grilles étendues 
   struct, border, indice, malloc, free, at(), swap. */
typedef struct ext_grid_t {
    double *grid_read,
	*left_read,
        *top_read,
        *bottom_read,
        *right_read,
        *leftR,
        *topR,
        *bottomR,
        *rightR,
	*grid_write,
	*left_write,
	*top_write,
	*bottom_write,
	*right_write;
    unsigned m_Nx,
	m_Ny;
} ext_grid_t;


enum border_t { TOP=1, BOTTOM=2, LEFT=4, RIGHT=8, NONE=0 };

enum request_t { TOP_REQUEST=0, BOTTOM_REQUEST=1, 
                 LEFT_REQUEST=2, RIGHT_REQUEST=3 };

unsigned IND(unsigned nx, unsigned ny, unsigned Ny) {
    return nx*Ny+ny;
}

/* fonction pour allouer automatiquement la mémoire à une grille étendue */
bool malloc_ext_grid(ext_grid_t* m_grid, bool with_write, bool with_borders) {
    bool alloc_state;
    m_grid->grid_read   = (double*) malloc((m_grid->m_Nx - 2) * (m_grid->m_Ny - 2) * sizeof(double));
    m_grid->top_read    = (double*) malloc(m_grid->m_Ny * sizeof (double));
    m_grid->bottom_read = (double*) malloc(m_grid->m_Ny * sizeof (double));  
    m_grid->left_read   = (double*) malloc(m_grid->m_Nx * sizeof (double));
    m_grid->right_read  = (double*) malloc(m_grid->m_Nx * sizeof (double));
    alloc_state = m_grid->grid_read && m_grid->top_read &&
	m_grid->bottom_read && m_grid->left_read && m_grid->right_read;	
    if (with_write) {
	m_grid->grid_write   = (double*) malloc((m_grid->m_Nx - 2) * (m_grid->m_Ny - 2) * sizeof(double));
	m_grid->top_write    = (double*) malloc(m_grid->m_Ny * sizeof (double));
	m_grid->bottom_write = (double*) malloc(m_grid->m_Ny * sizeof (double));  
	m_grid->left_write   = (double*) malloc(m_grid->m_Nx * sizeof (double));
	m_grid->right_write  = (double*) malloc(m_grid->m_Nx * sizeof (double));
	alloc_state = alloc_state && m_grid->grid_write && m_grid->top_write &&
	    m_grid->bottom_write && m_grid->left_write && m_grid->right_write;	
    }else{
        m_grid->grid_write   = NULL;
	m_grid->top_write    = NULL;
	m_grid->bottom_write = NULL;  
	m_grid->left_write   = NULL;
	m_grid->right_write  = NULL;
    }
    if (with_borders) {
	m_grid->topR    = (double*) malloc(m_grid->m_Ny * sizeof (double));
	m_grid->bottomR = (double*) malloc(m_grid->m_Ny * sizeof (double));
	m_grid->leftR   = (double*) malloc(m_grid->m_Nx * sizeof (double));
	m_grid->rightR  = (double*) malloc(m_grid->m_Nx * sizeof (double));
	alloc_state = alloc_state && m_grid->topR &&
	    m_grid->bottomR && m_grid->leftR && m_grid->rightR;	
    }else{
        m_grid->topR    = NULL;
	m_grid->bottomR = NULL;
	m_grid->leftR   = NULL;
	m_grid->rightR  = NULL;
    }
    return !alloc_state;
}


void free_ext_grid(ext_grid_t* ext_grid) {
    // PAS DE FREE DE EXT_GRID (pas alloué dynamiquement)
    
    if (ext_grid->grid_read)    free(ext_grid->grid_read);
    if (ext_grid->top_read)     free(ext_grid->top_read);
    if (ext_grid->bottom_read)  free(ext_grid->bottom_read);
    if (ext_grid->left_read)    free(ext_grid->left_read);
    if (ext_grid->right_read)   free(ext_grid->right_read);
    if (ext_grid->grid_write)   free(ext_grid->grid_write);
    if (ext_grid->top_write)    free(ext_grid->top_write);
    if (ext_grid->bottom_write) free(ext_grid->bottom_write);
    if (ext_grid->left_write)   free(ext_grid->left_write);
    if (ext_grid->right_write)  free(ext_grid->right_write);
    if (ext_grid->topR)         free(ext_grid->topR);
    if (ext_grid->bottomR)      free(ext_grid->bottomR);
    if (ext_grid->leftR)        free(ext_grid->leftR);
    if (ext_grid->rightR)       free(ext_grid->rightR);
}

/* lire facilement une grille étendue */
double* grid_read_at(const ext_grid_t* ext_grid, unsigned nx, unsigned ny) {

    // On a envoyé nx-1
    if (nx==-1) return &(ext_grid->topR[ny]);
    
    // On a envoyé nx+1
    if (nx==ext_grid->m_Nx) return &(ext_grid->bottomR[ny]);
    
    // On a envoyé nx-1 ou nx
    if (nx==0 && ny>=0 && ny<ext_grid->m_Ny) {
        return &(ext_grid->top_read[ny]);
    }
    // On a envoyé nx+1 ou nx    
    if (nx==ext_grid->m_Nx-1 && ny>=0 && ny<ext_grid->m_Ny) {
        return &(ext_grid->bottom_read[ny]);
    }
    if (ny==-1)return &(ext_grid->leftR[nx]); 
    if (ny==0) return &(ext_grid->left_read[nx]);
    if (ny==ext_grid->m_Ny-1) return &(ext_grid->right_read[nx]);
    if (ny==ext_grid->m_Ny) return &(ext_grid->rightR[nx]);
    return &(ext_grid->grid_read[(nx-1)*(ext_grid->m_Ny-2)+(ny-1)]);
}

/* écrire facilement dans une grille étendue */
void grid_write_at(const ext_grid_t* ext_grid, unsigned nx, unsigned ny, double val) {
    if (nx>0 && nx<ext_grid->m_Nx-1 && ny>0 && ny<ext_grid->m_Ny-1)
	ext_grid->grid_write[(nx-1)*(ext_grid->m_Ny-2)+(ny-1)] = val;
    else {
	if (nx==0) ext_grid->top_write[ny] = val;
	else ext_grid->bottom_write[ny] = val; // nx forcément égal à m_Nx-1
	if (ny==0) ext_grid->left_write[nx] = val;
	else ext_grid->right_write[nx] = val; // ny frocément égal = m_Ny-1
    }
}

/* échanger deux pointeurs de double */
void swap_dpts(double** a, double** b) {
    double* tmp=*a;
    *a = *b;
    *b = tmp;
}

/* échanger les parties read et write d'une grille étendue */
void swap_read_write(ext_grid_t* ext_grid) {
    swap_dpts(&(ext_grid->grid_read),   &(ext_grid->grid_write));
    swap_dpts(&(ext_grid->top_read),    &(ext_grid->top_write));
    swap_dpts(&(ext_grid->bottom_read), &(ext_grid->bottom_write));
    swap_dpts(&(ext_grid->left_read),   &(ext_grid->left_write));
    swap_dpts(&(ext_grid->right_read),  &(ext_grid->right_write));
}


/* Equations physiques et autres fonctions pour la résolution */
double isoline(double x, double y) {
    double theta = 3.141592653589793/2; 
    return sin(theta)*x+cos(theta)*y;
}

double exact_sol(double t, double x, double y,
		double T_max, double kappa, double sigma) {
    const double x_t = isoline(x, y);
    return T_max/sqrt( 1 + 4*t*kappa/pow(sigma, 2) ) *
	exp( -pow(x_t, 2)/(pow(sigma, 2) + 4*t*kappa) );
}

/* remplie notre grille locale à partir de ses coordonnées globales */
void exact_sol_grid(ext_grid_t* m_grid,
		    unsigned offset_nx, unsigned offset_ny, double dx, double dy,
		    double t, double T_max, double kappa, double sigma) {
    unsigned nx, ny;
    for (nx=0; nx<m_grid->m_Nx; ++nx)
	for (ny=0; ny<m_grid->m_Ny; ++ny) 
	    *grid_read_at(m_grid, nx, ny) = exact_sol(t, (offset_nx+nx)*dx,
						      (offset_ny+ny)*dy, T_max,
						      kappa, sigma);
}

/* résolution sur un point */  
void one_point(ext_grid_t* m_grid, unsigned border,
	       unsigned nx, unsigned ny,
	       double dt, double dx, double dy, double kappa) {
    if((border & TOP) && nx == 0) {
	m_grid->top_write[ny] = m_grid->top_read[ny];
	return;
    }
    if((border & BOTTOM) && nx == m_grid->m_Nx-1) {
	m_grid->bottom_write[ny] = m_grid->bottom_read[ny];
	return;
    }
    if((border & LEFT) && ny == 0) {
	m_grid->left_write[nx] = m_grid->left_read[nx];
	return;
    }
    if((border & RIGHT) && ny == m_grid->m_Ny-1) {
	m_grid->right_write[nx] = m_grid->right_read[nx];
	return;
    }

    double top=0, bottom=0, left=0, right=0, self=0;
    top    = *grid_read_at(m_grid, nx-1, ny);
    bottom = *grid_read_at(m_grid, nx+1, ny);
    left   = *grid_read_at(m_grid, nx, ny-1);
    right  = *grid_read_at(m_grid, nx, ny+1);
    self   = *grid_read_at(m_grid, nx, ny);
    
    grid_write_at(m_grid, nx, ny, self + dt * kappa *
    		  ( (left+right-2*self)/pow(dx, 2) + (top+bottom-2*self)/pow(dy, 2) ));
    return;
}

/* faire tousl les points 1 par 1 */
void one_step(ext_grid_t* m_grid, unsigned border,
	      double dt, double dx, double dy, double kappa) {
    unsigned nx, ny;
    for (nx=0; nx<m_grid->m_Nx; ++nx) 
	for (ny=0; ny<m_grid->m_Ny; ++ny) 
	    one_point(m_grid, border, 
		      nx, ny, dt, dx, dy,
		      kappa); 
}

/* diff L2 entre deux fonctions de R2 */
double L2_diff(const ext_grid_t* ext_grid_a, const ext_grid_t* ext_grid_b) {
    double diff = 0;
    unsigned nx, ny;
    if (ext_grid_a->m_Nx!=ext_grid_b->m_Nx || ext_grid_a->m_Ny!=ext_grid_b->m_Ny)
	return -1;
    for (ny=0; ny<ext_grid_a->m_Ny; ++ny) 
	for (nx=0; nx<ext_grid_a->m_Nx; ++nx) {
	    diff += pow((*grid_read_at(ext_grid_a, nx, ny)) -
			(*grid_read_at(ext_grid_b, nx, ny)), 2);
	}
    return diff;
}

/* Heuristique pour optimiser le nombre de threads utilisables 
   ne bug pas mais on pourrait faire mieux (EXPLE: pour 21 thr on 
   en a 1 unused mais si on prend 20 thr, les calculs changent et 
   au lieu de tous les utiliser on en a 2 unused*/
void optimize_Nx_thr_Ny_thr(int* Nx_threads_nb, int* Ny_threads_nb, int N_threads) {
    int diff[5];
    diff[0] = N_threads - (*Nx_threads_nb)*(*Ny_threads_nb);
    diff[1] = N_threads - (*Nx_threads_nb+1)*(*Ny_threads_nb);
    diff[2] = N_threads - (*Nx_threads_nb)*(*Ny_threads_nb+1);
    diff[3] = N_threads - (*Nx_threads_nb+1)*(*Ny_threads_nb-1);
    diff[4] = N_threads - (*Nx_threads_nb-1)*(*Ny_threads_nb+1);
    int min_ind=0,
	min=diff[0],
	ind;
    for (ind=1; ind<5; ind++) {
	if (diff[ind]<min && diff[ind]>=0) {
	    min_ind = ind;
	    min = diff[ind];
	}
    }
    if (min_ind==1) (*Nx_threads_nb)++;
    else if (min_ind==2) (*Ny_threads_nb)++;
    else if (min_ind==3) {
	(*Nx_threads_nb)++;
	(*Ny_threads_nb)--;
    }
    else if (min_ind==4) {
	(*Nx_threads_nb)--;
	(*Ny_threads_nb)++;
    }
}

void init_ext_grid_size(int rank, int N_threads, int Nx, int Ny,
			int* Nx_threads_nb, int* Ny_threads_nb, int* m_Nx, int* m_Ny,
			int* ny_thr, int* nx_thr,
                        int* offset_nx, int* offset_ny, unsigned* m_border) {

    // Les petits rectangles sont des "homotéties" du grand 
    *Nx_threads_nb = floor(sqrt(N_threads * Nx/Ny));
    *Ny_threads_nb = floor(sqrt(N_threads * Ny/Nx));
    // Heuristique sur la dimension des grilles étendues
    optimize_Nx_thr_Ny_thr(Nx_threads_nb, Ny_threads_nb, N_threads); 
    
    // On charge plus la majorité des threads pour qu'une minorité ait 
    // moins de travail (il ne faut pas non plus trop les surcharger).
    if(!(Nx % (*Nx_threads_nb - 1)) && !(Ny % (*Ny_threads_nb - 1))){
        *m_Nx = floor(Nx/(*Nx_threads_nb - 1));
        *m_Ny = floor(Ny/(*Ny_threads_nb - 1));
    }else{
        *m_Nx = floor(Nx/(*Nx_threads_nb));
        *m_Ny = floor(Ny/(*Ny_threads_nb));
    } 
    	   
    // Position de notre bout de grille sur la grille réelle
    *nx_thr = (rank/(*Ny_threads_nb));
    *ny_thr = (rank%(*Ny_threads_nb));
    *offset_nx = (*nx_thr) * (*m_Nx);
    *offset_ny = (*ny_thr) * (*m_Ny);

    if ((*nx_thr)==(*Nx_threads_nb) || (*ny_thr)==(*Ny_threads_nb)) {
	*m_Nx=0; *m_Ny=0;
    }
    //Si l'on se trouve sur le bord droit, on prend toutes les cases de droites restantes
    if((rank + 1) % (*Ny_threads_nb) == 0)
        *m_Ny = Ny - (*m_Ny) * (*Ny_threads_nb - 1); 
    // Idem si l'on se retrouve sur le bord du bas, on prend toutes les cases du bas
    if(rank/(*Ny_threads_nb) + 1 == *Nx_threads_nb)
        *m_Nx = Nx - (*m_Nx) * (*Nx_threads_nb - 1);

    // Choix du type de frontière
    *m_border = NONE;
    if(*offset_nx == 0) *m_border |= TOP;
    else if(*offset_nx + *m_Nx == Nx) *m_border |= BOTTOM;
    if(*offset_ny == 0) *m_border |= LEFT;
    else if(*offset_ny + *m_Ny == Ny) *m_border |= RIGHT;
}

// On attend que les envois et réceptions se soient bien passés.

void mpi_wait_requests(const unsigned m_border, 
                       MPI_Request* send_requests, 
                       MPI_Request* recv_requests){
    
    if(!(m_border & RIGHT)) {
        MPI_Wait(&send_requests[RIGHT_REQUEST], MPI_STATUS_IGNORE);
        MPI_Wait(&recv_requests[RIGHT_REQUEST], MPI_STATUS_IGNORE);
    }
    if(!(m_border & LEFT)) {
        MPI_Wait(&send_requests[LEFT_REQUEST], MPI_STATUS_IGNORE);
        MPI_Wait(&recv_requests[LEFT_REQUEST], MPI_STATUS_IGNORE);
    }
    if(!(m_border & TOP)) {
        MPI_Wait(&send_requests[TOP_REQUEST], MPI_STATUS_IGNORE);
        MPI_Wait(&recv_requests[TOP_REQUEST], MPI_STATUS_IGNORE);        
    }
    if(!(m_border & BOTTOM)) {
        MPI_Wait(&send_requests[BOTTOM_REQUEST], MPI_STATUS_IGNORE);
        MPI_Wait(&recv_requests[BOTTOM_REQUEST], MPI_STATUS_IGNORE);        
    }
}


int main(int argc, char **argv) {
    /*
      lecture des arguments
    */
    if (argc<13) {
	printf("Missing args (got %d instead of 12) : \n", argc);
        printf("Args : Nx, Ny, nb_steps, t_end, a, b, kappa, CFL, T_max, sigma, plot_every, calc_err\n");
	return 1;
    }
    /* 
       Initialisation de MPI 
    */
    MPI_Init(&argc, &argv);
    int N_threads, m_rank, tag = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &N_threads);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
    MPI_Request send_requests[4], recv_requests[4];

    /* 
       Définition des constantes 
    */
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
      Initialisation des variables
    */
    double *grid_tmp = NULL, *L2_err = NULL, *times = NULL;
    double L2_err_loc = 0,
	L2_err_global = 0,
        exec_time     = 0;
    
    int Nx_threads_nb, Ny_threads_nb,
	nx_thr, ny_thr,         // position de thr sur la grille
	m_Nx, m_Ny,             // Nb de points de discrétisation selon Ox et Oy inclus dans notre grilles
        offset_nx, offset_ny;   // Position de notre grille étendue dans la grille
    unsigned int m_border; 

    init_ext_grid_size(m_rank, N_threads, Nx, Ny,
		       &Nx_threads_nb, &Ny_threads_nb,
		       &m_Nx, &m_Ny, &nx_thr, &ny_thr,
		       &offset_nx, &offset_ny, &m_border);
    /* 
       Decompte des threads morts et création de deux nouveaux MPI_Comm adaptés
    */
    int N_alive_threads = Nx_threads_nb*Ny_threads_nb;
    int N_aux_alive_threads = ((N_threads > N_alive_threads) ? N_alive_threads+1 : N_alive_threads);
    unsigned rank_active_aux_thread = (N_aux_alive_threads == N_alive_threads ?  0 : N_alive_threads);
    int *alive_threads = (int*) malloc((N_alive_threads)*sizeof(int)),
	*aux_alive_threads = (int*) malloc(N_aux_alive_threads*sizeof(int)),
	athr;
    MPI_Group mpi_group_alive_world,
              mpi_group_aux_alive_world,
              mpi_group_world;
    MPI_Comm mpi_comm_alive_world, mpi_comm_aux_alive_world;
    for (athr=0; athr<N_alive_threads; ++athr){
	alive_threads[athr] = athr;
        aux_alive_threads[athr] = athr;
    }
    if(rank_active_aux_thread != 0) aux_alive_threads[N_aux_alive_threads-1] = N_aux_alive_threads-1;
    MPI_Comm_group(MPI_COMM_WORLD, &mpi_group_aux_alive_world);
    MPI_Comm_group(MPI_COMM_WORLD, &mpi_group_world);
    MPI_Group_incl(mpi_group_aux_alive_world, N_aux_alive_threads, aux_alive_threads, &mpi_group_aux_alive_world);
    MPI_Group_incl(mpi_group_world, N_alive_threads, alive_threads, &mpi_group_alive_world);
    MPI_Comm_create(MPI_COMM_WORLD, mpi_group_aux_alive_world, &mpi_comm_aux_alive_world);
    MPI_Comm_create(MPI_COMM_WORLD, mpi_group_alive_world, &mpi_comm_alive_world);

    if ((!m_Nx || !m_Ny) && m_rank != rank_active_aux_thread) { //Quitter si rien à faire --> regarder si MPI_Finalize ne gêne pas les autres threads
	printf("#rank = %d -- Nothing to do\n", m_rank);
	MPI_Finalize();
	return 0;
    }
    if (m_rank == rank_active_aux_thread) 
	printf("#rank = %d -- Auxiliary thread\n", m_rank);
    if (!m_rank)
	printf("\n#m_Nx = %d \n#m_Ny = %d \n#Nx_threads_nb = %d \n#Ny_threads_nb = %d\n",
	       m_Nx, m_Ny, Nx_threads_nb, Ny_threads_nb);
    
    /* 
       Initialisation de la grille étendue locale 
    */
    // Allocation de mémoire pour le centre et les bords de notre grille
    ext_grid_t m_grid, grid_sol_exact;
    bool alloc_err = false;
    if (m_Nx && m_Ny) {
	m_grid.m_Nx = m_Nx; grid_sol_exact.m_Nx = m_Nx;
	m_grid.m_Ny = m_Ny; grid_sol_exact.m_Ny = m_Ny;
	alloc_err = (malloc_ext_grid(&m_grid, true, true) || malloc_ext_grid(&grid_sol_exact, false, false));
    }
    
    // Si des threads ne font rien, il vont recevoir les erreurs L2
    if (calc_err && m_rank == rank_active_aux_thread) {
        L2_err = (double*) malloc(nb_steps * sizeof (double));
        times  = (double*) malloc(nb_steps * sizeof (double));
        alloc_err |= !times || !L2_err;
    } 
    if (alloc_err){
        printf("Cannot allocate memory (rank = %d)\n", m_rank);
	MPI_Finalize();
        return 2;
    }
    const double dt = t_end/nb_steps,
	dx = a/Nx,
	dy = b/Ny,
	CFL = kappa*dt/MIN(pow(dx, 2), pow(dy, 2));
    
    if(m_rank == rank_active_aux_thread) printf("CFL : %f \n", CFL); 
    
    if (CFL > CFL_max) {
        if (m_rank == rank_active_aux_thread) {
	    printf("CFL too big, aborting! \n");
	    fflush(stdout);
	}
	MPI_Finalize();
        return 3;
    } 
    /*
      initialisation des grilles locales (chaque thread initialise sa grille étendue locale)
    */
    if(m_rank != rank_active_aux_thread || !m_rank){
        exact_sol_grid(&m_grid, offset_nx, offset_ny, dx, dy,
                       0, T_max, kappa, sigma);
        exact_sol_grid(&grid_sol_exact, offset_nx, offset_ny, dx, dy,
                       0, T_max, kappa, sigma);
    }
    
    /*
      boucle principale
    */
    unsigned step;
    for (step=1; step<=nb_steps; step++) {
        
        if(m_rank != rank_active_aux_thread || !m_rank) {

            // On calcule la solution en t = dt*nb_step.
            one_step(&m_grid, m_border, dt, dx, dy, kappa);
            // On met la nouvelle solution dans la grille de lecture.
            swap_read_write(&m_grid);
            
            // Envoi à droite
            if(!(m_border & RIGHT))
                MPI_Isend(m_grid.right_read, m_Nx, MPI_DOUBLE, m_rank+1, tag, mpi_comm_alive_world, &send_requests[RIGHT_REQUEST]);
            if(!(m_border & LEFT))
                MPI_Irecv(m_grid.leftR, m_Nx, MPI_DOUBLE, m_rank-1, tag, mpi_comm_alive_world, &recv_requests[LEFT_REQUEST]);
            //Envoi à gauche
            if(!(m_border & LEFT))
                MPI_Isend(m_grid.left_read, m_Nx, MPI_DOUBLE, m_rank-1, tag, mpi_comm_alive_world, &send_requests[LEFT_REQUEST]);
            if(!(m_border & RIGHT))
                MPI_Irecv(m_grid.rightR, m_Nx, MPI_DOUBLE, m_rank+1, tag, mpi_comm_alive_world, &recv_requests[RIGHT_REQUEST]);
            //Envoi en haut
            if(!(m_border & TOP))
                MPI_Isend(m_grid.top_read, m_Ny, MPI_DOUBLE, m_rank - Ny_threads_nb, tag, mpi_comm_alive_world, &send_requests[TOP_REQUEST]);
            if(!(m_border & BOTTOM))
                MPI_Irecv(m_grid.bottomR, m_Ny, MPI_DOUBLE, m_rank + Ny_threads_nb, tag, mpi_comm_alive_world, &recv_requests[BOTTOM_REQUEST]);
            //Envoi en bas
            if(!(m_border & BOTTOM))
                MPI_Isend(m_grid.bottom_read, m_Ny, MPI_DOUBLE, m_rank + Ny_threads_nb, tag, mpi_comm_alive_world, &send_requests[BOTTOM_REQUEST]);
            if(!(m_border & TOP))
                MPI_Irecv(m_grid.topR, m_Ny, MPI_DOUBLE, m_rank - Ny_threads_nb, tag, mpi_comm_alive_world, &recv_requests[TOP_REQUEST]);
        }
        
        if (calc_err) {
	    if (m_rank != rank_active_aux_thread || !m_rank) {
		exact_sol_grid(&grid_sol_exact, offset_nx, offset_ny, dx, dy,
			   step*dt, T_max, kappa, sigma);
		L2_err_loc = L2_diff(&m_grid, &grid_sol_exact);
	    }
	    MPI_Reduce(&L2_err_loc, &L2_err_global, 1, MPI_DOUBLE, MPI_SUM, rank_active_aux_thread, mpi_comm_aux_alive_world);
	    if (m_rank == rank_active_aux_thread) {
		L2_err[step] = sqrt(L2_err_global) / (Nx * Ny);
		times[step] = step*dt;
		if (!(step%100) || step == 1) {
		    printf("\n Time, Error : %g, %g \n", times[step], L2_err[step]);
		    fflush(stdout);
		}
	    }
	}

        // On attend que tous les threads aient reçu et envoyé
        // leurs bords.
        if(m_rank != rank_active_aux_thread || !m_rank) 
            mpi_wait_requests(m_border, send_requests, recv_requests);
    }
    MPI_Barrier(mpi_comm_aux_alive_world);
    if(m_rank==rank_active_aux_thread) {
        exec_time = ((double) clock()/CLOCKS_PER_SEC);
        printf("\n Temps d'exécution du programme : %f s\n", exec_time);
    }
    
    /* Affichage graphique de l'erreur en fonction du temps */
    if (calc_err && m_rank==rank_active_aux_thread) {
	FILE* gnuplot_pipe = popen("gnuplot -p", "w");
	plot_curve_1D(gnuplot_pipe, times, L2_err,
		      nb_steps, "", "");
	pclose(gnuplot_pipe);
    }
    
    /* Libération des ressources */
    if (m_Nx && m_Ny) {
     	free_ext_grid(&m_grid);
     	free_ext_grid(&grid_sol_exact);
    }
    if (calc_err && m_rank==rank_active_aux_thread) {
        free(times);
        free(L2_err);
    }
    MPI_Group_free(&mpi_group_world);
    MPI_Group_free(&mpi_group_alive_world);
    MPI_Group_free(&mpi_group_aux_alive_world);
    if(m_rank != rank_active_aux_thread || !m_rank) MPI_Comm_free(&mpi_comm_alive_world);
    MPI_Comm_free(&mpi_comm_aux_alive_world); 
	    
    MPI_Finalize();
    return 0;
}