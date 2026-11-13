/* Cuánto tarda una carrera. Incluye el simulador sin tocarlo (su main se
   renombra) y mide cada correr() con el reloj monotónico.
       gcc -std=c99 -O2 -Dmain=principal -o medir_tiempo medir_tiempo.c -lm */
#define _POSIX_C_SOURCE 199309L
#include <time.h>
#include "electraton.c"
#undef main
static double ahora(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec / 1e9; }
int main(void) {
    const double ks[] = { 1.0, 1.2 };
    for (int i = 0; i < 2; i++) for (int e = 0; e < 3; e++) {
        Param p = PARAM; p.peukert = ks[i];
        double mejor = 1e9, peor = 0;
        Resultado r;
        for (int k = 0; k < 5; k++) {
            double t0 = ahora(); r = correr(&p, e, NULL); double dt = ahora() - t0;
            if (dt < mejor) mejor = dt; if (dt > peor) peor = dt;
        }
        printf("Peukert %.1f %-10s %5.1f min de carrera: %.3f s a %.3f s (5 veces)\n", p.peukert, NOMBRE[e], r.tiempo / 60, mejor, peor);
    }
    return 0;
}
