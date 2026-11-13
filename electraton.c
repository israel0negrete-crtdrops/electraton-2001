/* El Electratón vuelve a correr.

   Una reimplementación experimental, en C99 y sin más que la biblioteca
   estándar, de los principios de control de un proyecto de 2001: el control
   de energía de un vehículo eléctrico de competencia, hecho entonces con un
   ATMEL 8535 en ensamblador. No es una traducción de aquel código, que no se
   conserva: es lo que dice el documento (arranque suave y un comparador de
   60 A), puesto contra un simulador.

   Todo lo físico es hipótesis de modelado: la masa, el motor, la batería de
   plomo-ácido, el efecto Peukert y el recorrido no se recuperaron del vehículo
   original. Están en PARAM, juntos, para que se vean y se puedan mover.

       cc -std=c99 -O2 -o electraton electraton.c -lm
       ./electraton              las tres estrategias, con Peukert 1.0 a 1.3
       ./electraton --csv DIR    además, las series para las gráficas
       ./electraton --prueba     las pruebas

   Determinista: no hay azar. Las tres estrategias corren el mismo recorrido,
   con el mismo piloto, la misma batería y los mismos parámetros. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------ hipótesis */

typedef struct {
    /* vehículo y piloto */
    double masa;          /* kg, con piloto */
    double crr;           /* resistencia a la rodadura */
    double cda;           /* m², área frontal por coeficiente de arrastre */
    double rueda;         /* m, radio */
    double reduccion;     /* vueltas del motor por vuelta de la rueda */
    double eta;           /* eficiencia de la transmisión */
    /* motor de corriente directa con imán permanente */
    double k;             /* V·s/rad = N·m/A */
    double ra;            /* ohm, armadura */
    double la;            /* H, inductancia de armadura */
    double rcable;        /* ohm, cables y transistor */
    /* batería: dos de 12 V de plomo-ácido en serie */
    double ah;            /* Ah a 20 horas */
    double peukert;       /* exponente */
    double rint;          /* ohm, las dos juntas */
    double voc_lleno;     /* V en vacío, cargada */
    double voc_vacio;     /* V en vacío, descargada */
    /* recorrido: dos rectas y dos curvas */
    double recta;         /* m */
    double curva;         /* m */
    double v_curva;       /* m/s, la que el piloto lleva en la curva */
    double frenado;       /* m/s², lo que frena el piloto */
    /* control */
    double rampa;         /* s, de 0 a 100 % de ciclo de trabajo */
    double limite;        /* A, el umbral del comparador */
    double paso_limite;   /* cuánto baja el ciclo de trabajo por milisegundo arriba del umbral */
    /* simulación */
    double dt;            /* s */
    double telemetria;    /* s, una muestra cada tanto, como en 2001 */
} Param;

static const Param PARAM = {
    .masa = 150, .crr = 0.010, .cda = 0.20, .rueda = 0.25, .reduccion = 4.0, .eta = 0.95,
    .k = 0.11, .ra = 0.06, .la = 0.0004, .rcable = 0.01,
    .ah = 35, .peukert = 1.2, .rint = 0.03, .voc_lleno = 25.5, .voc_vacio = 23.6,
    .recta = 350, .curva = 150, .v_curva = 7.0, .frenado = 2.5,
    .rampa = 3.0, .limite = 60, .paso_limite = 0.02,
    .dt = 0.001, .telemetria = 20,
};

#define G 9.81
#define RHO 1.2

typedef enum { DIRECTO = 0, RAMPA = 1, RAMPA_LIMITE = 2 } Estrategia;
static const char *NOMBRE[] = { "directo", "rampa", "rampa+60A" };

/* ------------------------------------------------------------ resultado */

typedef struct {
    double tiempo, distancia;          /* hasta vaciar la batería */
    double distancia30;                /* m a los 30 minutos: la carrera corta, donde la batería no se acaba */
    int vueltas;
    double vuelta_media;               /* s por vuelta, sin la primera */
    double pico_bat, pico_motor;       /* A, el mayor pulso (corriente del motor) de cada milisegundo */
    double ah, ah_efectivos, wh;       /* lo que salió de la batería */
    double perdidas_bat, perdidas_motor;   /* Wh en resistencias */
    double trabajo_motor;              /* Wh entregados al eje */
    double seg_sobre_limite;           /* s con la corriente del motor arriba del umbral */
    double pico_telemetria;            /* A, el mayor que ve una muestra cada 20 s */
    int muestras_sobre_limite;         /* muestras de 20 s arriba del umbral */
    double v_max;                      /* m/s */
    double v_min_curva;                /* m/s, la menor dentro de una curva, después de la primera vuelta */
    double balance;                    /* error relativo del balance de energía */
} Resultado;

/* ------------------------------------------------------------ el recorrido y el piloto */

/* Posición dentro de la vuelta: 0..recta es recta, luego curva, recta, curva. */
static int en_curva(const Param *p, double s, double *falta_curva) {
    double vuelta = 2 * (p->recta + p->curva);
    double x = fmod(s, vuelta);
    double tramo = p->recta + p->curva;
    double y = fmod(x, tramo);
    if (y < p->recta) { *falta_curva = p->recta - y; return 0; }
    *falta_curva = 0;
    return 1;
}

/* El piloto: en recta, acelerador a fondo hasta el punto de frenado de la
   curva; en la curva, mantiene la velocidad de la curva. El mismo piloto
   para las tres estrategias. Devuelve el acelerador (0..1) y el freno (N). */
static void piloto(const Param *p, double s, double v, double *acelerador, double *freno) {
    double falta;
    *freno = 0;
    if (!en_curva(p, s, &falta)) {
        if (v * v > p->v_curva * p->v_curva + 2 * p->frenado * falta) {
            *acelerador = 0;
            *freno = p->masa * p->frenado;
        } else *acelerador = 1;
    } else {
        double e = p->v_curva - v;
        *acelerador = e > 0 ? fmin(1, 3.0 * e) : 0;
        if (v > p->v_curva + 0.5) *freno = p->masa * p->frenado;
    }
}

/* ------------------------------------------------------------ una carrera */

typedef struct { FILE *serie, *muestras; double cada; } Salida;

static Resultado correr(const Param *p, Estrategia e, const Salida *out) {
    Resultado r;
    memset(&r, 0, sizeof r);
    double t = 0, s = 0, v = 0, im = 0, ib = 0, ciclo = 0, tope = 1;
    double ah_ef = 0, e_voc = 0, e_mec = 0;
    double i20 = p->ah / 20;
    double vuelta = 2 * (p->recta + p->curva);
    int vuelta_actual = 0;
    double t_vuelta1 = 0, t_ultima = 0;
    double sig_muestra = p->telemetria, sig_serie = 0;
    double rmot = p->ra + p->rcable;
    double tau = p->la / rmot, decae = exp(-p->dt / tau);

    while (ah_ef < p->ah && t < 6 * 3600) {
        double acel_piloto, freno;
        piloto(p, s, v, &acel_piloto, &freno);

        /* el control: qué ciclo de trabajo pide */
        double pedido = acel_piloto;
        if (e != DIRECTO) {
            double sube = p->dt / p->rampa;
            pedido = acel_piloto < ciclo ? acel_piloto : fmin(acel_piloto, ciclo + sube);
        }
        if (e == RAMPA_LIMITE) {
            /* el comparador: arriba del umbral, el tope baja; abajo, se recupera a ritmo de rampa */
            if (im > p->limite) tope = fmax(0, tope - p->paso_limite);
            else tope = fmin(1, tope + p->dt / p->rampa);
            pedido = fmin(pedido, tope);
        }
        ciclo = pedido;

        /* la batería, con la corriente del paso anterior. Sin capacitor que la alise: mientras el
           transistor conduce, la batería da toda la corriente del motor; mientras no, nada. Es la
           hipótesis que menos favorece al control, que pasa más tiempo con ciclo parcial. */
        double soc = 1 - ah_ef / p->ah;
        double voc = p->voc_vacio + (p->voc_lleno - p->voc_vacio) * soc;
        double vt = voc - p->rint * im;

        /* el motor: la corriente de armadura tiende a su valor de equilibrio con la constante L/R;
           con un solo transistor no puede invertirse */
        double w = v / p->rueda * p->reduccion;
        double iss = (ciclo * vt - p->k * w) / rmot;
        im = iss + (im - iss) * decae;
        if (im < 0) im = 0;
        ib = ciclo * im;                  /* el promedio de lo que sale de la batería */

        /* la mecánica */
        double f_motor = p->k * im * p->reduccion * p->eta / p->rueda;
        double f_res = p->masa * G * p->crr + 0.5 * RHO * p->cda * v * v;
        double f = f_motor - (v > 0 ? f_res + freno : 0);
        double a = f / p->masa;
        if (v <= 0 && a < 0) a = 0;
        double v_nueva = fmax(0, v + a * p->dt);
        s += 0.5 * (v + v_nueva) * p->dt;
        v = v_nueva;
        t += p->dt;

        /* la contabilidad */
        /* Peukert sobre los pulsos: una fracción `ciclo` del tiempo, a la corriente del motor */
        double ief = im > 0 ? ciclo * im * pow(im / i20, p->peukert - 1) : 0;
        ah_ef += ief * p->dt / 3600;
        r.ah += ib * p->dt / 3600;
        r.wh += vt * ib * p->dt / 3600;
        e_voc += voc * ib * p->dt;
        r.perdidas_bat += ciclo * p->rint * im * im * p->dt / 3600;
        r.perdidas_motor += rmot * im * im * p->dt / 3600;
        e_mec += p->k * im * w * p->dt;
        if (im > r.pico_bat) r.pico_bat = im;     /* el pulso: lo que da la batería mientras conduce */
        if (im > r.pico_motor) r.pico_motor = im;
        if (im > p->limite) r.seg_sobre_limite += p->dt;
        if (v > r.v_max) r.v_max = v;
        if (r.distancia30 == 0 && t >= 1800) r.distancia30 = s;
        { double f; if (s > vuelta && en_curva(p, s, &f) && (r.v_min_curva == 0 || v < r.v_min_curva)) r.v_min_curva = v; }

        int n = (int)(s / vuelta);
        if (n > vuelta_actual) { vuelta_actual = n; t_ultima = t; if (n == 1) t_vuelta1 = t; }

        if (t >= sig_muestra - 1e-9) {
            /* una muestra instantánea cada 20 s, como la telemetría de 2001 */
            if (im > r.pico_telemetria) r.pico_telemetria = im;
            if (im > p->limite) r.muestras_sobre_limite++;
            if (out && out->muestras) fprintf(out->muestras, "%s,%.0f,%.2f,%.2f,%.2f,%.3f\n", NOMBRE[e], t, im, ib, vt, v);
            sig_muestra += p->telemetria;
        }
        /* la serie fina, sólo los primeros tres minutos */
        if (out && out->serie && t <= 180 && t >= sig_serie - 1e-9) {
            fprintf(out->serie, "%s,%.3f,%.2f,%.2f,%.3f,%.3f,%.3f\n", NOMBRE[e], t, ib, im, ciclo, v, vt);
            sig_serie += out->cada;
        }
    }
    r.tiempo = t;
    r.distancia = s;
    r.vueltas = vuelta_actual;
    r.vuelta_media = vuelta_actual > 1 ? (t_ultima - t_vuelta1) / (vuelta_actual - 1) : 0;
    r.ah_efectivos = ah_ef;
    r.trabajo_motor = e_mec / 3600;
    /* lo que dio la fuerza electromotriz de la batería = eje + resistencias */
    double entra = e_voc / 3600, sale = r.trabajo_motor + r.perdidas_motor + r.perdidas_bat;
    r.balance = fabs(entra - sale) / entra;
    return r;
}

/* ------------------------------------------------------------ salida */

static void imprimir(const Param *p, Estrategia e, const Resultado *r) {
    printf("  %-10s  %5.1f min  %6.2f km (%5.2f km a 30 min)  %3d vueltas  %5.1f s/vuelta  pico %5.0f A (motor %5.0f A)  "
           "%5.0f s arriba de %.0f A  telemetría ve %5.1f A  %5.1f Ah (%5.1f efectivos)  %6.0f Wh  "
           "pérdidas bat %4.0f Wh motor %4.0f Wh  vmax %4.1f km/h\n",
           NOMBRE[e], r->tiempo / 60, r->distancia / 1000, r->distancia30 / 1000, r->vueltas, r->vuelta_media, r->pico_bat, r->pico_motor,
           r->seg_sobre_limite, p->limite, r->pico_telemetria, r->ah, r->ah_efectivos, r->wh,
           r->perdidas_bat, r->perdidas_motor, r->v_max * 3.6);
}

static void json(FILE *f, double k, Estrategia e, const Resultado *r, int ultimo) {
    fprintf(f, "  {\"peukert\": %.2f, \"estrategia\": \"%s\", \"minutos\": %.2f, \"km\": %.3f, \"km_30min\": %.3f, \"vueltas\": %d, "
               "\"s_por_vuelta\": %.2f, \"pico_bateria_A\": %.1f, \"pico_motor_A\": %.1f, \"s_sobre_60A\": %.1f, "
               "\"pico_telemetria_A\": %.1f, \"muestras_sobre_60A\": %d, \"ah\": %.2f, \"ah_efectivos\": %.2f, "
               "\"wh\": %.1f, \"perdidas_bateria_wh\": %.1f, \"perdidas_motor_wh\": %.1f, \"trabajo_motor_wh\": %.1f, "
               "\"vmax_kmh\": %.2f}%s\n",
            k, NOMBRE[e], r->tiempo / 60, r->distancia / 1000, r->distancia30 / 1000, r->vueltas, r->vuelta_media, r->pico_bat, r->pico_motor,
            r->seg_sobre_limite, r->pico_telemetria, r->muestras_sobre_limite, r->ah, r->ah_efectivos, r->wh,
            r->perdidas_bat, r->perdidas_motor, r->trabajo_motor, r->v_max * 3.6, ultimo ? "" : ",");
}

/* ------------------------------------------------------------ pruebas */

static int fallas = 0;
#define CHECA(c, ...) do { if (!(c)) { fallas++; printf("FALLA: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

static void pruebas(void) {
    Param p = PARAM;
    Resultado d = correr(&p, DIRECTO, NULL), ra = correr(&p, RAMPA, NULL), li = correr(&p, RAMPA_LIMITE, NULL);
    /* la energía cuadra: lo que da la batería es eje más resistencias */
    CHECA(d.balance < 0.01 && ra.balance < 0.01 && li.balance < 0.01, "balance %.4f %.4f %.4f", d.balance, ra.balance, li.balance);
    /* el comparador: la corriente del motor no pasa del umbral más que un poco */
    CHECA(li.pico_motor < p.limite * 1.15, "pico con límite %.1f", li.pico_motor);
    /* sin control, el arranque pide más de 60 A */
    CHECA(d.pico_bat > p.limite, "pico directo %.1f", d.pico_bat);
    /* con Peukert 1, los Ah efectivos son los Ah */
    { Param q = p; q.peukert = 1.0; Resultado x = correr(&q, DIRECTO, NULL);
      CHECA(fabs(x.ah - x.ah_efectivos) < 1e-6 * x.ah, "peukert 1: %.4f %.4f", x.ah, x.ah_efectivos); }
    /* determinista */
    { Resultado x = correr(&p, RAMPA, NULL); CHECA(x.distancia == ra.distancia, "determinista"); }
    /* la telemetría cada 20 s nunca ve más que el pico real */
    CHECA(d.pico_telemetria <= d.pico_bat && li.pico_telemetria <= li.pico_bat, "telemetría");
    /* el piloto sostiene la curva: dentro de ella no baja de 6.5 m/s ni frena de más, con ninguna estrategia */
    CHECA(d.v_min_curva > 6.5 && ra.v_min_curva > 6.5 && li.v_min_curva > 6.5, "curva %.2f %.2f %.2f", d.v_min_curva, ra.v_min_curva, li.v_min_curva);
    /* las cifras de la nota (Peukert 1.2), redondeadas como en la nota */
    #define R1(x) (round((x) * 10) / 10)
    CHECA(R1(d.distancia / 1000) == 27.6 && R1(ra.distancia / 1000) == 28.3 && R1(li.distancia / 1000) == 30.6, "km hasta vaciar");
    CHECA(R1(d.distancia30 / 1000) == 17.4 && R1(ra.distancia30 / 1000) == 17.3 && R1(li.distancia30 / 1000) == 16.9, "km en 30 min");
    CHECA(round(d.pico_bat) == 254 && round(ra.pico_bat) == 191 && round(li.pico_bat) == 64, "picos");
    CHECA(round(d.seg_sobre_limite) == 258 && d.muestras_sobre_limite == 12 && round(d.pico_telemetria) == 110, "telemetría de la nota");
    /* sin pérdidas de cobertura ni azar: la vuelta media es la misma en una segunda corrida */
    { Resultado x = correr(&p, DIRECTO, NULL); CHECA(x.vuelta_media == d.vuelta_media, "vuelta media"); }
    printf(fallas ? "%d fallas\n" : "ok\n", fallas);
}

/* ------------------------------------------------------------ todo */

int main(int argc, char **argv) {
    if (argc > 1 && !strcmp(argv[1], "--prueba")) { pruebas(); return fallas != 0; }
    const char *dir = (argc > 2 && !strcmp(argv[1], "--csv")) ? argv[2] : NULL;
    const double ks[] = { 1.0, 1.1, 1.2, 1.3 };
    FILE *fj = NULL;
    Salida out = { NULL, NULL, 0.01 };
    char ruta[1024];
    if (dir) {
        snprintf(ruta, sizeof ruta, "%s/resultados.json", dir); fj = fopen(ruta, "w");
        snprintf(ruta, sizeof ruta, "%s/serie_3min.csv", dir); out.serie = fopen(ruta, "w");
        snprintf(ruta, sizeof ruta, "%s/telemetria_20s.csv", dir); out.muestras = fopen(ruta, "w");
        if (!fj || !out.serie || !out.muestras) { perror(dir); return 1; }
        fprintf(out.serie, "estrategia,t_s,i_bateria_promedio_A,i_pulso_A,ciclo,v_ms,v_bateria_V\n");
        fprintf(out.muestras, "estrategia,t_s,i_pulso_A,i_bateria_promedio_A,v_bateria_V,v_ms\n");
        fprintf(fj, "[\n");
    }
    for (int i = 0; i < 4; i++) {
        Param p = PARAM;
        p.peukert = ks[i];
        printf("Peukert %.1f\n", p.peukert);
        for (int e = 0; e < 3; e++) {
            /* las series sólo con la hipótesis de inicio */
            Resultado r = correr(&p, e, (dir && ks[i] == PARAM.peukert) ? &out : NULL);
            imprimir(&p, e, &r);
            if (fj) json(fj, p.peukert, e, &r, i == 3 && e == 2);
        }
    }
    /* el compromiso: el umbral del comparador, de 30 A a 200 A, con la hipótesis de inicio */
    printf("Barrido del umbral (Peukert %.1f, rampa + límite)\n", PARAM.peukert);
    FILE *fb = NULL;
    if (dir) { snprintf(ruta, sizeof ruta, "%s/barrido_umbral.csv", dir); fb = fopen(ruta, "w");
               if (fb) fprintf(fb, "umbral_A,km_hasta_vaciar,minutos,km_30min,s_por_vuelta,pico_motor_A\n"); }
    const double umbrales[] = { 30, 40, 50, 60, 80, 100, 150, 200 };
    for (int i = 0; i < 8; i++) {
        Param p = PARAM;
        p.limite = umbrales[i];
        Resultado r = correr(&p, RAMPA_LIMITE, NULL);
        printf("  %3.0f A: %6.2f km en %5.1f min, %5.2f km a 30 min, %5.1f s/vuelta, pico motor %5.0f A\n",
               p.limite, r.distancia / 1000, r.tiempo / 60, r.distancia30 / 1000, r.vuelta_media, r.pico_motor);
        if (fb) fprintf(fb, "%.0f,%.3f,%.2f,%.3f,%.2f,%.1f\n", p.limite, r.distancia / 1000, r.tiempo / 60, r.distancia30 / 1000, r.vuelta_media, r.pico_motor);
    }
    if (fb) fclose(fb);
    if (fj) { fprintf(fj, "]\n"); fclose(fj); fclose(out.serie); fclose(out.muestras); }
    return 0;
}
