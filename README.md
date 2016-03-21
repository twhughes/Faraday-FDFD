# Faraday-FDFD

This is an implementation of the finite difference frequency domain algorithm (FDFD) for electromagnetics written in Julia


The solver and it's helper functions live in 'src' with the main field solver being fdfd.jl


This function can be called as follows

        (Ex,Ey,Ez,Hx,Hy,Hz) = fdfd(ER2,MUR2,RES,NPML,BC,lambda0,Pol,Q; verbose=false,TFSF=false,theta=0);


        User inputs:
                ER2        array                        relative permittivity distribution on 2x grid               complex array
                MUR2       array                        relative permeability distribution on 2x grid               complex array
                RES        [dx dy]                      grid resolution / size of dx, dy on 2x grid                 positive float array
                NPML       [Nxlo Nxhi Nylo Nyhi]        pml border sizes (in grid points)                           non-negative integer array
                lambda0    float                        free space incident wavelength                              positive float
                Pol        {'Ez','Hz'}                  polarization (Ez,Hx,Hy) or (Hz,Ex,Ey)                       string
                Q          array                        if TFSF, the masking matrix for SF. else, source Jz or Hz   complex array
        Options:
                verbose    {true,false}                 print out timings and status
                TFSF       {true,false}                 use Q array as TFSF constructor (true) or as source (false)
                theta      float                        angle of incidence (degrees) of TFSF source                 float {-90 .. 90}

Due to plotting issues with Julia, it works best to do post / pre processing on a Jupyter notebook running julia.
Included are some example notebooks.
I recommend SimpleCalc.ipynb to get started.


This is a work in progress.. come back later for added functionality
