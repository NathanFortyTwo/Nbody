/*
** nbody_brute_force.c - nbody simulation using the brute-force algorithm (O(n*n))
** Cuda version
**/

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>
#include <unistd.h>
#include <cuda.h> 
#include <cuda_runtime.h>

#ifdef DISPLAY
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#endif

#include "ui.h"
#include "nbody.h"
#include "cuda_stuff.cuh"

FILE* f_out = NULL;

int nparticles = 10; /* number of particles */
float T_FINAL = 1.0; /* simulation end time */
particle_t* particles;

double sum_speed_sq = 0;
double max_acc = 0;
double max_speed = 0;

/** CUDA **/
cudaError_t err;
particle_t* d_particles;
dim3 threadsPerBlock;
dim3 numBlocks;

void init()
{
  /** Initialisation des variables nbthreadbyblock et nbblockbygrid **/
  threadsPerBlock = dim3(256);
  numBlocks = dim3((nparticles + threadsPerBlock.x - 1) / threadsPerBlock.x);
  printf("Threads per block: %d\n", threadsPerBlock.x);
  printf("Blocks per grid: %d\n", numBlocks.x);
}


#ifdef DISPLAY
extern Display* theDisplay; /* These three variables are required to open the */
extern GC theGC;            /* particle plotting window.  They are externally */
extern Window theMain;      /* declared in ui.h but are also required here.   */
#endif

/* compute the force that a particle with position (x_pos, y_pos) and mass 'mass'
 * applies to particle p
 * Runs as a kernel on the GPU
 */
__device__ __host__ void compute_force(particle_t* p, double x_pos, double y_pos, double mass) {
  double x_sep, y_sep, dist_sq, grav_base;

  x_sep = x_pos - p->x_pos;
  y_sep = y_pos - p->y_pos;
  dist_sq = MAX((x_sep * x_sep) + (y_sep * y_sep), 0.01);

  /* Use the 2-dimensional gravity rule: F = d * (GMm/d^2) */
  grav_base = GRAV_CONSTANT * (p->mass) * (mass) / dist_sq;

  p->x_force += grav_base * x_sep;
  p->y_force += grav_base * y_sep;
}

/* compute all forces that all particles apply to all others
 * Runs as a kernel on the GPU
 */
__global__ void compute_forces_kernel(particle_t* d_particles, int nparticles)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < nparticles)
  {
    int j;
    d_particles[i].x_force = 0;
    d_particles[i].y_force = 0;
    for (j = 0; j < nparticles; j++)
    {
      particle_t* p = &d_particles[j];
      /* compute the force of particle j on particle i */
      compute_force(&d_particles[i], p->x_pos, p->y_pos, p->mass);
    }
  }
}

/* compute the new position/velocity */
void move_particle(particle_t* p, double step) {

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
void all_move_particles(double step)
{
  cudaMemcpy(d_particles, particles, sizeof(particle_t) * nparticles, cudaMemcpyHostToDevice);
  /** Lancement du kernel **/
  compute_forces_kernel << <numBlocks, threadsPerBlock >> > (d_particles, nparticles);
  cudaMemcpy(particles, d_particles, sizeof(particle_t) * nparticles, cudaMemcpyDeviceToHost);

  /* then move all particles and return statistics */
  int i;
  for (i = 0; i < nparticles; i++)
  {
    move_particle(&particles[i], step);
  }
}


void print_all_particles(FILE* f)
{
  int i;
  for (i = 0; i < nparticles; i++)
  {
    particle_t* p = &particles[i];
    fprintf(f, "particle={pos=(%f,%f), vel=(%f,%f)}\n", p->x_pos, p->y_pos, p->x_vel, p->y_vel);
  }
}

void run_simulation()
{
  double t = 0.0, dt = 0.01;
  while (t < T_FINAL && nparticles > 0)
  {
    /* Update time. */
    t += dt;
    printf("t = %lf\n", t);
    /* Move particles with the current and compute rms velocity. */
    all_move_particles(dt);

    /* Adjust dt based on maximum speed and acceleration--this
       simple rule tries to insure that no velocity will change
       by more than 10% */

    dt = 0.1 * (max_speed) / (max_acc);

  }
}

/*
  Place particles in their initial positions.
*/
void all_init_particles(int num_particles, particle_t* particles)
{
  int    i;
  double total_particle = num_particles;

  for (i = 0; i < num_particles; i++) {
    particle_t* particle = &particles[i];

    particle->x_pos = i * 2.0 / nparticles - 1.0;
    particle->y_pos = 0.0;
    particle->x_vel = 0.0;
    particle->y_vel = particle->x_pos;

    particle->mass = 1.0 + (num_particles + i) / total_particle;
    particle->node = NULL;

  }
}

/*
  Simulate the movement of nparticles particles.
*/
int main(int argc, char** argv)
{
  if (argc >= 2)
  {
    nparticles = atoi(argv[1]);
  }
  if (argc == 3)
  {
    T_FINAL = atof(argv[2]);
  }

  init();

  /* Allocate global shared arrays for the particles data set. */
  /** Allocation memoire sur le host(CPU) **/
  particles = (particle_t*)malloc(sizeof(particle_t) * nparticles);
  all_init_particles(nparticles, particles);

  /** Allocation memoire sur le device(GPU) **/
  err = cudaMalloc((void**)&d_particles, sizeof(particle_t) * nparticles);
  gpuErrchk(err);

  /** Transfert mémoire du host vers le device **/
  err = cudaMemcpy(d_particles, particles, sizeof(particle_t) * nparticles, cudaMemcpyHostToDevice);
  gpuErrchk(err);

  /* Initialize thread data structures */
#ifdef DISPLAY
  /* Open an X window to display the particles */
  simple_init(100, 100, DISPLAY_SIZE, DISPLAY_SIZE);
#endif

  struct timeval t1, t2;
  gettimeofday(&t1, NULL);

  run_simulation();

  gettimeofday(&t2, NULL);

  double duration = (t2.tv_sec - t1.tv_sec) + ((t2.tv_usec - t1.tv_usec) / 1e6);


#ifdef DUMP_RESULT
  FILE* f_out = fopen("particles.log", "w");
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
  /** Libération de la mémoire **/
  err = cudaFree(d_particles);
  gpuErrchk(err);
  free(particles);

  return 0;
}
