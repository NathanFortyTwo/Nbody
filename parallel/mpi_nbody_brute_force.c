/*
** nbody_brute_force.c - nbody simulation using the brute-force algorithm (O(n*n))
**
**/

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>
#include <unistd.h>

#ifdef DISPLAY
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#endif

#include "ui.h"
#include "nbody.h"
#include "nbody_tools.h"
#include <mpi.h>
#include <omp.h>
FILE *f_out = NULL;

int nparticles = 10; /* number of particles */
float T_FINAL = 1.0; /* simulation end time */
particle_t *particles;

double sum_speed_sq = 0;
double max_acc = 0;
double max_speed = 0;

void init()
{
  /* Nothing to do */
}

#ifdef DISPLAY
extern Display *theDisplay; /* These three variables are required to open the */
extern GC theGC;            /* particle plotting window.  They are externally */
extern Window theMain;      /* declared in ui.h but are also required here.   */
#endif

/* compute the force that a particle with position (x_pos, y_pos) and mass 'mass'
 * applies to particle p
 */
void compute_force(particle_t *p, double x_pos, double y_pos, double mass)
{
  double x_sep, y_sep, dist_sq, grav_base;

  x_sep = x_pos - p->x_pos;
  y_sep = y_pos - p->y_pos;
  dist_sq = MAX((x_sep * x_sep) + (y_sep * y_sep), 0.01);

  /* Use the 2-dimensional gravity rule: F = d * (GMm/d^2) */
  grav_base = GRAV_CONSTANT * (p->mass) * (mass) / dist_sq;

  p->x_force += grav_base * x_sep;
  p->y_force += grav_base * y_sep;
}

/* compute the new position/velocity */
void move_particle(particle_t *p, double step)
{

  p->x_pos += (p->x_vel) * step;
  p->y_pos += (p->y_vel) * step;
  double x_acc = p->x_force / p->mass;
  double y_acc = p->y_force / p->mass;
  p->x_vel += x_acc * step;
  p->y_vel += y_acc * step;

  /* compute statistics */
  double cur_acc = (x_acc * x_acc + y_acc * y_acc);
  cur_acc = sqrt(cur_acc);
  double speed_sq = (p->x_vel) * (p->x_vel) + (p->y_vel) * (p->y_vel);
  double cur_speed = sqrt(speed_sq);

  sum_speed_sq += speed_sq;
  max_acc = MAX(max_acc, cur_acc);
  max_speed = MAX(max_speed, cur_speed);
}

/*
  Move particles one time step.

  Update positions, velocity, and acceleration.
  Return local computations.
*/
void all_move_particles(double step, int rank, int size)
{
  /* First calculate force for particles. */
  int i;
  // #pragma omp parallel for shared(particles)

  int start = rank * (nparticles / size);
  int end = (rank + 1) * (nparticles / size);

  // if (rank == size - 1)
  //{
  //   end = nparticles;
  // }

  // on suppose pour l'instant que nparticles est divisible par size
  assert(end - start == nparticles / size);

  for (i = start; i < end; i++)
  {
    int j;
    particles[i].x_force = 0;
    particles[i].y_force = 0;
    for (j = 0; j < nparticles; j++)
    {
      particle_t *p = &particles[j];
      /* compute the force of particle j on particle i */
      compute_force(&particles[i], p->x_pos, p->y_pos, p->mass);
    }
  }
  double xforces[end - start];
  double yforces[end - start];
  for (i = start; i < end; i++)
  {
    xforces[i - start] = particles[i].x_force;
    yforces[i - start] = particles[i].y_force;
  }

  // gather all forces
  MPI_Barrier(MPI_COMM_WORLD);

  int displs[size];
  for (i = 0; i < size; i++)
  {
    displs[i] = i * (nparticles / size);
  }

  double allxforces[nparticles];
  double allyforces[nparticles];

  // int recvcounts[size];
  // for (i = 0; i < size; i++)
  //{
  //   recvcounts[i] = end - start;
  // }

  // MPI_Allgatherv(xforces, end - start, MPI_DOUBLE, allxforces, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
  // MPI_Allgatherv(yforces, end - start, MPI_DOUBLE, allyforces, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

  MPI_Allgather(xforces, end - start, MPI_DOUBLE, allxforces, end - start, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(yforces, end - start, MPI_DOUBLE, allyforces, end - start, MPI_DOUBLE, MPI_COMM_WORLD);

  for (i = start; i < end; i++)
  {
    particles[i].x_force = allxforces[i - start];
    particles[i].y_force = allyforces[i - start];
  }

  /* then move all particles and return statistics */
  // #pragma omp parallel for shared(particles)
  for (i = 0; i < nparticles; i++)
  {
    move_particle(&particles[i], step);
  }
}

/* display all the particles */
void draw_all_particles()
{
  int i;
  for (i = 0; i < nparticles; i++)
  {
    int x = POS_TO_SCREEN(particles[i].x_pos);
    int y = POS_TO_SCREEN(particles[i].y_pos);
    draw_point(x, y);
  }
}

void print_all_particles(FILE *f)
{
  int i;
  for (i = 0; i < nparticles; i++)
  {
    particle_t *p = &particles[i];
    fprintf(f, "particle={pos=(%f,%f), vel=(%f,%f)}\n", p->x_pos, p->y_pos, p->x_vel, p->y_vel);
  }
}

void run_simulation(int rank, int size)
{
  double t = 0.0, dt = 0.01;
  while (t < T_FINAL && nparticles > 0)
  {
    /* Update time. */
    t += dt;
    /* Move particles with the current and compute rms velocity. */
    all_move_particles(dt, rank, size);

    /* Adjust dt based on maximum speed and acceleration--this
       simple rule tries to insure that no velocity will change
       by more than 10% */

    dt = 0.1 * max_speed / max_acc;

    /* Plot the movement of the particle */
#if DISPLAY
    clear_display();
    draw_all_particles();
    flush_display();
#endif
  }
}

/*
  Simulate the movement of nparticles particles.
*/
int main(int argc, char **argv)
{
  if (argc >= 2)
  {
    nparticles = atoi(argv[1]);
  }
  if (argc == 3)
  {
    T_FINAL = atof(argv[2]);
  }

  // MPI_Init(&argc, &argv);
  // int rank;
  // MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  init();
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /* Allocate global shared arrays for the particles data set. */
  particles = malloc(sizeof(particle_t) * nparticles);
  all_init_particles(nparticles, particles);

  /* Initialize thread data structures */
#ifdef DISPLAY
  /* Open an X window to display the particles */
  simple_init(100, 100, DISPLAY_SIZE, DISPLAY_SIZE);
#endif

  struct timeval t1, t2;
  gettimeofday(&t1, NULL);

  /* Main thread starts simulation ... */
  run_simulation(rank, size);
  if (rank)
  {
    return 0;
  }
  gettimeofday(&t2, NULL);

  double duration = (t2.tv_sec - t1.tv_sec) + ((t2.tv_usec - t1.tv_usec) / 1e6);

#ifdef DUMP_RESULT
  FILE *f_out = fopen("particles.log", "w");
  assert(f_out);
  print_all_particles(f_out);
  fclose(f_out);
#endif

  printf("-----------------------------\n");
  printf("nparticles: %d\n", nparticles);
  printf("T_FINAL: %f\n", T_FINAL);
  printf("-----------------------------\n");
  printf("Simulation took %lf s to complete\n", duration);

#ifdef DISPLAY
  clear_display();
  draw_all_particles();
  flush_display();

  printf("Hit return to close the window.");

  getchar();
  /* Close the X window used to display the particles */
  XCloseDisplay(theDisplay);
#endif
  // MPI_Finalize();
  return 0;
}
