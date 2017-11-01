#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <stdbool.h>
#include <mpi.h>
#define MIN(a, b) ((a)<(b) ? (a) : (b))

/* Définition et fonctions sur les grilles étendues 
   struct, border, indice, malloc, free, at(), swap. */
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
    unsigned m_Nx,
	m_Ny;
} ext_grid_t;

enum border_t { TOP=1, BOTTOM=2, LEFT=4, RIGHT=8, NONE=0 };

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
    }
    if (with_borders) {
	m_grid->topR    = (double*) malloc(m_grid->m_Ny * sizeof (double));
	m_grid->bottomR = (double*) malloc(m_grid->m_Ny * sizeof (double));
	m_grid->leftR   = (double*) malloc(m_grid->m_Nx * sizeof (double));
	m_grid->rightR  = (double*) malloc(m_grid->m_Nx * sizeof (double));
	alloc_state = alloc_state && m_grid->topR &&
	    m_grid->bottomR && m_grid->leftR && m_grid->rightR;	
    }
    return !alloc_state;
}

/* fonction pour free la mémoire -> NE PAS UTILISER, cause des segfaults,
   peut être à cause d'autres erreurs... ;) */
void free_ext_grid(ext_grid_t* ext_grid) {
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
    free(ext_grid);
}

/* lire facilement une grille étendue -> le bug est probablement la dedans... 
   Il faut vérifier la dernière ligne tout particulièrement. Même si ça pas 
   l'air d'être celle la qui plante. Le printf qui pollue la sortie est 
   juste en dessous la vvvvv */
double* grid_read_at(const ext_grid_t* ext_grid, unsigned nx, unsigned ny) {
    printf("%d %d\n", nx, ny); fflush(stdout);
    if (nx==-1) return &(ext_grid->topR[ny]);
    if (nx==0) return &(ext_grid->top_read[ny]);
    if (nx==ext_grid->m_Nx-1) return &(ext_grid->bottom_read[ny]);
    if (nx==ext_grid->m_Nx) return &(ext_grid->bottomR[ny]);

    if (ny==-1)return &(ext_grid->leftR[nx]); 
    if (ny==0) return &(ext_grid->left_read[nx]);
    if (ny==ext_grid->m_Ny-1) return &(ext_grid->right_read[nx]);
    if (ny==ext_grid->m_Ny) return &(ext_grid->rightR[nx]);
    return &(ext_grid->grid_read[(nx-1)*(ext_grid->m_Ny-2)+(ny-1)]);
}

/* écrire facilement dans une grille étendue -> le bug est peut être aussi ici... */
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

/* échanger deux pointeurs, c'est assez simple et ça à l'air de marcher... 
   Cf code commenté après la fonction main pour vérifier */
void swap_dpts(double** a, double** b) {
    double* tmp=*a;
    *a = *b;
    *b = tmp;
}

/* échanger les parties read et write d'une grille étendue, ça à l'air de marcher aussi */
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

/* remplie notre grille locale à partir de ses coordonnées globales -> semble OK */
void exact_sol_grid(ext_grid_t* m_grid,
		    unsigned offset_nx, unsigned offset_ny, double dx, double dy,
		    double t, double T_max, double kappa, double sigma) {
    unsigned nx, ny;
    for (ny=0; ny<m_grid->m_Ny; ++ny) 
	for (nx=0; nx<m_grid->m_Ny; ++nx) 
	    *grid_read_at(m_grid, nx, ny) = exact_sol(t, (offset_nx+nx)*dx,
						      (offset_ny+ny)*dy, T_max,
						      kappa, sigma);
}

/* résolution sur un point -> d'après mes test il arrive pas à lire les 
   coord (nx, ny) = (m_Nx-1, -1) càd left pour le thread rank=1. Cependant
   si tu le mets tout seul avant genre printf grid_read_at(&m_grid, 369, -1) 
   il marche sans problème... C'est la que je me suis arrêté... Et il 
   n'y a pas de problème en rank=0 c'est parce que border=LEFT donc ça passe 
   pas par la fonction grid_read_at */  
void one_point(ext_grid_t* m_grid, enum border_t border,
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

    //printf("t"); fflush(stdout);
    top    = *grid_read_at(m_grid, nx-1, ny);
    //printf("b"); fflush(stdout);
    bottom = *grid_read_at(m_grid, nx+1, ny);
    //printf("l"); fflush(stdout);
    left   = *grid_read_at(m_grid, nx, ny-1);
    //printf("r"); fflush(stdout);
    right  = *grid_read_at(m_grid, nx, ny+1);
    //printf("s"); fflush(stdout);
    self   = *grid_read_at(m_grid, nx, ny);
    //printf("w\n"); fflush(stdout);
    grid_write_at(m_grid, nx, ny, self + dt * kappa *
    		  ( (left+right-2*self)/pow(dx, 2) + (top+bottom-2*self)/pow(dy, 2) ));
    return;
}

/* faire tousl les points 1 par 1 */
void one_step(ext_grid_t* m_grid, enum border_t border,
	      double dt, double dx, double dy, double kappa) {
    unsigned nx, ny;
    for (nx=0; nx<m_grid->m_Nx; ++nx) {
	for (ny=0; ny<m_grid->m_Ny; ++ny) {
	    //printf("%d %d %d\n", nx, ny, border); fflush(stdout);
	    one_point(m_grid, border, 
		      nx, ny, dt, dx, dy,
		      kappa); 
	}
    }
}

/* diff L2 entre deux fonctions de R2 */
double L2_diff(const ext_grid_t* ext_grid_a, const ext_grid_t* ext_grid_b) {
    double diff = 0;
    unsigned nx, ny;
    if (ext_grid_a->m_Nx!=ext_grid_b->m_Nx || ext_grid_a->m_Ny!=ext_grid_b->m_Ny)
	return -1;
    for (ny=0; ny<ext_grid_a->m_Ny; ++ny) 
	for (nx=0; nx<ext_grid_a->m_Nx; ++nx) {
	    if ( *grid_read_at(ext_grid_a, nx, ny)>1) printf("%d %d\n", nx, ny, *grid_read_at(ext_grid_a, nx, ny), *grid_read_at(ext_grid_b, nx, ny));
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
    for (ind=1; ind<5; ind++) 
	if (diff[ind]<min && diff[ind]>0) {
	    min_ind = ind;
	    min = diff[ind];
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
                        int* offset_nx, int* offset_ny, enum border_t* m_border) {

    // Les petits rectangles sont des "homotéties" du grand 
    *Nx_threads_nb = floor(sqrt(N_threads * Nx/Ny));
    *Ny_threads_nb = floor(sqrt(N_threads * Ny/Nx));
    // Heuristique sur la dimension des grilles étendues
    optimize_Nx_thr_Ny_thr(Nx_threads_nb, Ny_threads_nb, N_threads); 
    
    *m_Nx = floor(Nx/(*Nx_threads_nb));
    *m_Ny = floor(Ny/(*Ny_threads_nb));
	   
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

int main(int argc, char **argv) {
    /*
      lecture des arguments
    */
    if (argc<13) {
	printf("Missing args (got %d instead of 12) : \n", argc);
        printf("Args : Nx, Ny, nb_steps, t_end, a, b, kappa, CFL, T_max, sigma, plot_every, calc_err\n");
	return EXIT_FAILURE;
    }
    /* 
       Initialisation de MPI 
    */
    MPI_Init(&argc, &argv);
    int N_threads, m_rank, tag = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &N_threads);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
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
	L2_err_global = 0;
    
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
       Decompte des threads morts et création d'un nouveau MPI_Comm adapté
    */
    int N_alive_threads = Nx_threads_nb*Ny_threads_nb;
    int * alive_threads = (int*) malloc((N_alive_threads)*sizeof(int));
    int athr;
    MPI_Group mpi_group_alive_world,
	mpi_group_world;
    MPI_Comm mpi_comm_alive_world;
    for (athr=0; athr<N_alive_threads; ++athr)
	alive_threads[athr] = athr;
    MPI_Comm_group(MPI_COMM_WORLD, &mpi_group_world);
    MPI_Group_incl(mpi_group_world, N_alive_threads, alive_threads, &mpi_group_alive_world);
    MPI_Comm_create(MPI_COMM_WORLD, mpi_group_alive_world, &mpi_comm_alive_world);

    if (!m_Nx || !m_Ny) { //Quitter si rien à faire --> regarder si MPI_Finalize ne gêne pas les autres threads
	printf("#rank = %d -- Nothing to do\n", m_rank);
	MPI_Finalize();
	return EXIT_SUCCESS;
    }
    if (!m_rank)
	printf("\n#m_Nx = %d \n#m_Ny = %d \n#Nx_threads_nb = %d \n#Ny_threads_nb = %d \n#Unused threads = %d\n",
	       m_Nx, m_Ny, Nx_threads_nb, Ny_threads_nb, N_threads-(Nx_threads_nb*Ny_threads_nb));

    /* 
       Initialisation de la grille étendue locale 
    */
    ext_grid_t m_grid, grid_sol_exact;
    m_grid.m_Nx = m_Nx; grid_sol_exact.m_Nx = m_Nx;
    m_grid.m_Ny = m_Ny; grid_sol_exact.m_Ny = m_Ny;
    // Allocation de mémoire pour le centre et les bords de notre grille
    bool alloc_err = (malloc_ext_grid(&m_grid, true, true) || malloc_ext_grid(&grid_sol_exact, false, false));

    // Le thread 0 reçoit les erreurs des autres threads et calcule l'erreur.
    if (calc_err && !m_rank) {
	L2_err = (double*) malloc(nb_steps*sizeof(double));
	times = (double*) malloc(nb_steps*sizeof(double));
        alloc_err |= !times || !L2_err;
    } 
    if (alloc_err){
        printf("Cannot allocate memory\n");
	MPI_Finalize();
        return EXIT_FAILURE;
    }	
    const double dt = t_end/nb_steps,
	dx = a/Nx,
	dy = b/Ny,
	CFL = kappa*dt/MIN(pow(dx, 2), pow(dy, 2));
    /* if(m_rank == 0) printf("CFL : %f \n", CFL); */
    /* if (CFL>CFL_max) { */
    /* 	if(m_rank == 0) printf("CFL too big... Aborting! \n"); */
    /* 	return EXIT_FAILURE; */
    /* } */
    /*
      initialisation des grilles locales (chaque thread initialise sa grille étendue locale)
      premier calcul d'erreur -> évident mais permet de débug deux trois trucs. Marche bien maintenant.
    */
    exact_sol_grid(&m_grid, offset_nx, offset_ny, dx, dy,
		   0, T_max, kappa, sigma);
    exact_sol_grid(&grid_sol_exact, offset_nx, offset_ny, dx, dy,
		   0, T_max, kappa, sigma);
    if(calc_err) {
	if(m_rank) {
	    L2_err_loc = L2_diff(&m_grid, &grid_sol_exact);
	    printf("%d %f\n", m_rank, L2_err_loc);
	    MPI_Reduce(&L2_err_loc, &L2_err_global, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_alive_world);
	} else {
	    L2_err_loc = L2_diff(&m_grid, &grid_sol_exact);
	    printf("%d %f\n", m_rank, L2_err_loc);
	    MPI_Reduce(&L2_err_loc, &L2_err_global, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_alive_world);
	    L2_err[0] = sqrt(L2_err_global)/(Nx*Ny);
	    times[0] = 0;
	    printf("Time : %f \n Error : %f \n", times[0], L2_err[0]);
	}
    }
    /*
      boucle principale
    */
    MPI_Barrier(mpi_comm_alive_world);
    unsigned step;
    for (step=1; step<nb_steps; step++) {
	/* 
	   Les comms marchent apparemment bien, càd si il n'y a qu'elles dans la boucle ça termine
	   Idée : Enlever one_step qui est compliquée et simplement essayer de lire intégralement 
	   chacune des grille étendue pour voir si on les lit bien déjà ou si il y pas des problèmes
	   d'allocation mémoire allouée mais caca quand même
	*/
        // Envoi à droite
        if(!(m_border & RIGHT))
            MPI_Send(m_grid.right_read, m_Nx, MPI_DOUBLE, m_rank+1, tag, mpi_comm_alive_world);
        if(!(m_border & LEFT))
            MPI_Recv(m_grid.leftR, m_Nx, MPI_DOUBLE, m_rank-1, tag, mpi_comm_alive_world, MPI_STATUS_IGNORE);
        //Envoi à gauche
	if(!(m_border & LEFT))
            MPI_Send(m_grid.left_read, m_Nx, MPI_DOUBLE, m_rank-1, tag, mpi_comm_alive_world);
        if(!(m_border & RIGHT))
            MPI_Recv(m_grid.rightR, m_Nx, MPI_DOUBLE, m_rank+1, tag, mpi_comm_alive_world, MPI_STATUS_IGNORE);
        //Envoi en haut
        if(!(m_border & TOP))
            MPI_Send(m_grid.top_read, m_Ny, MPI_DOUBLE, m_rank - Ny_threads_nb, tag, mpi_comm_alive_world);
        if(!(m_border & BOTTOM))
            MPI_Recv(m_grid.bottomR, m_Ny, MPI_DOUBLE, m_rank + Ny_threads_nb, tag, mpi_comm_alive_world, MPI_STATUS_IGNORE);
        //Envoi en bas
        if(!(m_border & BOTTOM))
            MPI_Send(m_grid.bottom_read, m_Ny, MPI_DOUBLE, m_rank + Ny_threads_nb, tag, mpi_comm_alive_world);
        if(!(m_border & TOP))
            MPI_Recv(m_grid.topR, m_Ny, MPI_DOUBLE, m_rank - Ny_threads_nb, tag, mpi_comm_alive_world, MPI_STATUS_IGNORE);
        // On avance d'un pas dans le temps

	MPI_Barrier(mpi_comm_alive_world);
	usleep(m_rank*5e6);
	printf("\n\nrang: %d -- pos: %d %d -- border %d\n", m_rank, nx_thr, ny_thr, m_border); fflush(stdout);
	printf("grid_read    %p\n", m_grid.grid_read);
	printf("top_read     %p\n", m_grid.top_read);
	printf("bottom_read  %p\n", m_grid.bottom_read);
	printf("left_read    %p\n", m_grid.left_read);
	printf("right_read   %p\n", m_grid.right_read);
	printf("grid_write   %p\n", m_grid.grid_write);
	printf("top_write    %p\n", m_grid.top_write);
	printf("bottom_write %p\n", m_grid.bottom_write);
	printf("left_write   %p\n", m_grid.left_write);
	printf("right_write  %p\n", m_grid.right_write);
	printf("topR         %p\n", m_grid.topR);
	printf("bottomR      %p\n", m_grid.bottomR);
	printf("leftR        %p\n", m_grid.leftR);
	printf("rightR       %p\n", m_grid.rightR);
	
	//if (m_rank==1) one_step(&m_grid, m_border, dt, dx, dy, kappa);
	swap_read_write(&m_grid);


        // Réception de toutes les erreurs sur le thread 0.
        if(calc_err){
             exact_sol_grid(&grid_sol_exact, offset_nx, offset_ny, dx, dy,
			    step*dt, T_max, kappa, sigma);
             if(m_rank) {
                 L2_err_loc = L2_diff(&m_grid, &grid_sol_exact);
                 MPI_Reduce(&L2_err_loc, &L2_err_global, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_alive_world);
             } else {
                 L2_err_loc = L2_diff(&m_grid, &grid_sol_exact);
		 printf("%d %f\n", m_rank, L2_err_loc);
                 MPI_Reduce(&L2_err_loc, &L2_err_global, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_alive_world);
                 L2_err[step] = sqrt(L2_err_global)/(Nx*Ny);
                 times[step] = step*dt;
                 printf("Time : %f \n Error : %f \n", times[step], L2_err[step]);
             }
        }
    }
    MPI_Finalize();
    return EXIT_SUCCESS;

}
       
    /* printf("\n#rank = %d \n#m_border = %u \n#m_Nx = %d \n#m_Ny = %d \n#Nx_threads_nb = %d \n#Ny_threads_nb = %d \n#nx_thr = %d \n#ny_thr = %d \n#offset_nx = %d \n#offset_ny = %d \n", */
    /* 	   m_rank, m_border, m_Nx, m_Ny, Nx_threads_nb, Ny_threads_nb, nx_thr, ny_thr, offset_nx, offset_ny); */
    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* unsigned nx, ny; */
    /* for (nx=0; nx<m_Nx; ++nx) */
    /* 	for (ny=0; ny<m_Ny; ++ny) { */
    /* 	    printf("%d %d\n", offset_nx+nx, offset_ny+ny); */
    /* 	    fflush(stdout); */
    /* 	} */

    /* usleep(m_rank*1e4); */
    /* printf("\nrank - %d\n", m_rank); */
    /* printf("grid_read    %p\n", m_grid.grid_read);    */
    /* printf("top_read     %p\n", m_grid.top_read);     */
    /* printf("bottom_read  %p\n", m_grid.bottom_read);  */
    /* printf("left_read    %p\n", m_grid.left_read);    */
    /* printf("right_read   %p\n", m_grid.right_read);   */
    /* printf("grid_write   %p\n", m_grid.grid_write);   */
    /* printf("top_write    %p\n", m_grid.top_write);    */
    /* printf("bottom_write %p\n", m_grid.bottom_write); */
    /* printf("left_write   %p\n", m_grid.left_write);   */
    /* printf("right_write  %p\n", m_grid.right_write);  */
    /* printf("topR         %p\n", m_grid.topR);         */
    /* printf("bottomR      %p\n", m_grid.bottomR);      */
    /* printf("leftR        %p\n", m_grid.leftR);        */
    /* printf("rightR       %p\n", m_grid.rightR);  */    
    
    /* free_ext_grid(&m_grid); */
    /* free_ext_grid(&grid_sol_exact); */
    /* if (calc_err) { */
    /* 	free(times); */
    /* 	free(L2_err); */
    /* } */
