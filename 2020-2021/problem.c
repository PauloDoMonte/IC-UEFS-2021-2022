#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "tools.h"

/*
	As unidades dessa simulação são kilometros, kilogramas e anos
*/

float soma = 0;
float q = 0;

double particulas[3][8] = {
//massa,raio,semieixo,ecentricidade,inclinacao,argumento do perigeu,longitude do no ascendente
  {1.989e+30,696340*1000}, // Sol
  {5.973332e+24,6378.1366,149.60e6,0.01671022,0.00005,-11.26064,102.94719,0.0},
  {1000,0.4,1.698665726802664*1.498e8,0.5077360606964316,12.16734013772536,224.6831644667536,356.6558403536075,0}
};

void heartbeat(struct reb_simulation* r){

  struct reb_particle p = r->particles[2];
  soma += sqrt((p.x*p.x)+(p.y*p.y)+(p.z*p.z))/1.498e8;
  q += 1;
  printf("Media: %f\n", soma/q);
  printf("T:%f\n", r->t);
  printf("%f\n", r->dt);

//  for (int i = 0; i < r->N; i++){
//    struct reb_particle p = r->particles[i];
//    printf("Particula%i\tDistancia:%f\n",i,(sqrt((p.x*p.x)+(p.y*p.y)+(p.z*p.z)))/1.498e8);
//  }
//  printf("%f\n", r->dt);
}

int main(int argc, char* argv[]){

  struct reb_simulation* r = reb_create_simulation();
  r->G=6.646596924499661e-05;
  r->integrator = REB_INTEGRATOR_IAS15;
  //r->dt=-1e-3;
  //r->ri_whfast.corrector=17;
  //r->ri_whfast.safe_mode=0;
  //r->ri_whfast.kernel = REB_WHFAST_KERNEL_LAZY;
  r->heartbeat=heartbeat;

  struct reb_particle primary = {0};
  primary.m = particulas[0][0];
  primary.r = particulas[0][1];
  reb_add(r,primary);
  
  struct reb_particle terra = reb_tools_orbit_to_particle(r->G, primary,particulas[1][0],particulas[1][2],particulas[1][3],particulas[1][4],particulas[1][6],particulas[1][5],0);
  reb_add(r,terra);

  struct reb_particle neo = reb_tools_orbit_to_particle(r->G, primary,particulas[2][0],particulas[2][2],particulas[2][3],particulas[2][4],particulas[2][6],particulas[2][5],0);
  reb_add(r,neo);

  reb_move_to_com(r);

  reb_integrate(r,-100000);

}

