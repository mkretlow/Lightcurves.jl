# Lightcurves.jl
Julia Package for Asteroid (Rotational) Light Curves and Phase Curves.


## Observed Lightcurves
Analyse rotational light curves of small bodies (asteroids).

Einbauen:

(G)LS,FALC,PDM,ANOVA

* https://www.cbabelgium.com/peranso/UserGuideHTML/Lomb-ScargleGLS.html
* https://juliaastro.github.io/LombScargle.jl/stable/#Introduction
* https://github.com/ckm3/Py-PDM

* Direkter Daten download von ALCDEF, sowie Ephemeride von MPC (für Delta,r,PHA usw.) oder mit 2-body-ephem berechnen.

* Plots für MPB
* Ausgabe der Daten im ALCDEF Format für upload

## Synthetic Lightcurves
Create synthetic LCs for different models and/or for small body 3D mesh files.


## Phasecurves
Tools to work on phasecurves of small bodies (asteroids, TNOs, etc.).

- min_pha([dat1,dat2]) : predict minimum solar phase angle for a list of small bodies

* Auswertung von beob. Phasenkurven => H,slope


## LC Inversion
* Inversion von Kaasaleinen implementieren (wrapper zum C code)
* Julia Interface to ADAM (tbd)

## Plots
* Measurement plots
* half split plot


# MEMO

https://stackoverflow.com/questions/34361250/fitting-time-series-with-fourier-components-estimating-fourier-series-coefficie

https://www.cbabelgium.com/peranso/index.html
