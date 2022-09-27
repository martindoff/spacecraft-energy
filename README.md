# spacecraft-energy
Spacecraft energy management

Source code for the spacecraft energy management algorithm presented in the paper: ["Spacecraft Energy Management Using Convex Optimisation"](https://www.researchgate.net/profile/Martin-Doff-Sotta/publication/352863815_Spacecraft_energy_management_using_convex_optimisation/links/60dcf81d92851ca9449b4a1f/Spacecraft-energy-management-using-convex-optimisation.pdf) 

# Gettting started

In the MATLAB command window, run

```
simulate
```

This will run the open loop optimisation.

If you want to run the closed-loop (MPC) optimisation instead, uncomment line 66 of `simulate.m`:


```
cvx_sat_MPC;
```

and comment line 67:

```
%ccvx_sat_open_loop;
```
