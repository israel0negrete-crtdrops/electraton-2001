# The Electratón runs again

An experimental reimplementation, in C99 with nothing but the standard library,
of the principles of a 2001 energy controller for an electric competition
vehicle (an *Electratón*). The original ran on an ATMEL 8535 in assembly; that
code is not kept, so this is **not a translation**. It takes what the 2001
document records—a soft start and a 60 A comparator—and puts it against a
simulator.

Note: [The Electratón runs again](https://crtdrops.mx/notas/the-electraton-runs-again.html)
(Spanish: [El Electratón vuelve a correr](https://crtdrops.mx/notas/el-electraton-vuelve-a-correr.html)).

## Everything physical is an assumption

The mass, the motor, the lead-acid battery, the Peukert exponent and the track
were **not** recovered from the original vehicle. They are all in `PARAM`, at
the top of `electraton.c`, so they can be seen and changed:

| | Assumption |
|---|---|
| Vehicle with driver | 150 kg, 0.25 m wheel, 4:1 reduction |
| Motor | permanent-magnet DC, k = 0.11 V·s/rad, 0.06 Ω, 0.4 mH |
| Battery | two 12 V lead-acid in series, 35 Ah (20 h), 0.03 Ω, Peukert 1.2 (also run at 1.0, 1.1, 1.3) |
| Track | 1 km lap: two 350 m straights, two 150 m curves at 7 m/s (about 25 km/h) |
| Driver | flat out on the straight, brakes for the curve, holds the curve speed; the same for every strategy |

## Strategies

1. **direct:** duty cycle = throttle.
2. **ramp:** the duty cycle rises at most 0→100 % in 3 s; it drops immediately.
3. **ramp+60A:** also, while the motor current is above 60 A, the duty cycle is pulled down.

Same track, same driver, same initial charge, same parameters. No randomness.

## Run

    cc -std=c99 -O2 -o electraton electraton.c -lm
    ./electraton --prueba          # tests
    ./electraton                   # summary, Peukert 1.0–1.3, and a sweep of the threshold
    ./electraton --csv resultados  # also writes resultados.json, serie_3min.csv,
                                   # telemetria_20s.csv and barrido_umbral.csv

Names and comments are in Spanish, as written.

## Results

`resultados/pi5/` holds the run on the lab's Raspberry Pi 5, with the
machine, the compiler and the timing in `equipo.txt`.

## What it does not model

The "operating point after braking" from 2001 (its definition was not
recovered); regenerative braking; battery temperature and recovery between
pulses; the motor's magnetic saturation; a driver who adapts to the
controller; the competition's real rules.

Peukert's law and the battery's resistive losses are applied to each current
pulse of the chopper (the battery gives the full motor current while the
transistor conducts, and nothing otherwise), with no bus capacitor smoothing
it. That is the assumption that least favors the current limit, which spends
more time at partial duty; with a capacitor, the limit's advantage would be
larger.

MIT license.
