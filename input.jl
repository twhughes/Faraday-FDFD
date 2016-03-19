###################################################################################################################################################
#
#   FDFD algorithm input
#   Written by Tyler Hughes
#   twhughes@stanford.edu
#   PhD Candidate, Stanford Dept. of Applied Physics.
#
#   User inputs:
#       ER2        array                        relative permittivity distribution on 2x grid       complex array
#       MUR2       array                        relative permeability distribution on 2x grid       complex array (same size as ER2)
#       RES        [dx dy]                      grid resolution / size of dx, dy on 2x grid         positive float array
#       NPML       [Nxlo Nxhi Nylo Nyhi]        pml border sizes (in grid points)                   non-negative integer array
#       lambda0    float                        free space incident wavelength                      positive float
#       Pol        {'Ez','Hz'}                  polarization (Ez,Hx,Hy) or (Hz,Ex,Ey)               string
#       theta      float                        angle of incidence (degrees)                        float {-90 .. 90}
#       Q          array                        scattered field mask array (1 where SF, 0 else)     boolean array (same size as ER2,MUR2)
###################################################################################################################################################


Nx = 20; Ny = 20;
ER2  = ones(Complex64,Nx,Ny);
MUR2 = ones(Complex64,Nx,Ny);
RES  = [1e-6,1e-6];
NPML = [0 0 0 0]
lambda0 = 1e-6;
Pol = 'Hz';
theta = 0;
Q = ones(Bool,Nx,Ny);;

[Ex,Ey,Ez,Hx,Hy,Hz] = fdfd(ER2,MUR2,RES,NPML,lambda0,Pol,theta,Q);
