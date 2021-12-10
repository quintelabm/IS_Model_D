# IS_Model_D

Immune System Model considering the Damage caused due to inflammation.

Continuation of the work published by [Quintela et al 2014](https://doi.org/10.1155/2014/410457)

## Mathematical Model

The mathematical model represents the dynamics of the activation of the adaptive immune response at the nearest lymph node via migration of an antigen presenting cell. The tissue is discretized as an hexahedron and the equations are solved with finite differences. To represent the flux between the two compartments it is considered that only the cells that are in contact with an area that is set to be a blood vessel or lymph vessel can migrate.

### Baseline Model: 

The variables considered to diffuse in the tissue are:
- Resting antigen presenting cells (RM)
- Activated antigen presenting cells (AM)
- Antigen (A)
- Antibodies (F)

The variables with dynamics considered within the lymph node:
- Antigen presenting cells that migrated (AM_L)
- T helper cells (Th)
- B cells (B)
- Plasma cells (P)
- Antibodies (F_L)

![model](https://static-01.hindawi.com/articles/bmri/volume-2014/410457/figures/410457.fig.003.jpg)


Representation of the blood vessels:

![bloodvessel](https://static-01.hindawi.com/articles/bmri/volume-2014/410457/figures/410457.fig.007b.jpg)

The lymph vessels were positioned differently:

![lymphvessels](https://static-01.hindawi.com/articles/bmri/volume-2014/410457/figures/410457.fig.0010b.jpg)

In addition to the model showed on the figure above, aka the baseline model, pro and anti-inflammatory cytokines and a new variable to represent the damage to the tissue were added.



### How to Run

```
./make
```


### Developed by 

Prof. Barbara de Melo Quintela

