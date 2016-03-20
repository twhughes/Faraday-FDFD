###################################################################################################################################################
#
#   FDFD algorithm.
#   Written by Tyler Hughes
#   twhughes@stanford.edu
#   PhD Candidate, Stanford Dept. of Applied Physics.
#
#   User inputs:
#       ER2        array                        relative permittivity distribution on 2x grid       complex array
#       MUR2       array                        relative permeability distribution on 2x grid       complex array
#       RES        [dx dy]                      grid resolution / size of dx, dy on 2x grid         positive float array
#       NPML       [Nxlo Nxhi Nylo Nyhi]        pml border sizes (in grid points)                   non-negative integer array
#       lambda0    float                        free space incident wavelength                      positive float
#       Pol        {'Ez','Hz'}                  polarization (Ez,Hx,Hy) or (Hz,Ex,Ey)               string
#       theta      float                        angle of incidence (degrees)                        float {-90 .. 90}
#
###################################################################################################################################################

include("calcpml_2d.jl");

function fdfd(ER2,MUR2,RES,NPML,BC,lambda0,Pol,theta,Q; verbose=false);

    if (verbose)
        tic();
        tic();
        println("(0) loading variables and doing error checking");
    end

    # (0) load variables and do error checking
    k0 = 2*pi/lambda0;
    (dx2,dy2) = RES;
    dx = dx2*2;  dy = dy2*2;
    (Nxlo,Nxhi,Nylo,Nyhi) = NPML;
    (Nx2,Ny2) = size(ER2);
    NGRID = (Nx2,Ny2);

    if (verbose)
        toc();
        tic();
        println("(1) determining material properties in reflected and transmitted regions");
    end
    # (1) determine material properties in reflected and transmitted regions

    e_ref = 1;   e_trans = 1;
    mu_ref = 1;  mu_trans = 1;
    n_ref = sqrt(mu_ref*e_ref);
    n_trans = sqrt(mu_trans*e_trans);

    if (verbose)
        toc();
        tic();
        println("(2) calculating PML terms");
    end

    # (2) compute PML terms s_x(x,y) s_y(x,y)
    (sx,sy) = calcpml_2d(NGRID,NPML);

    if (verbose)
        toc();
        tic();
        println("(3) incorporating PML into 2x material grid");
    end

    # (3) incorporate PML into 2x material grid
    ERxx = ER2./sx.*sy;
    ERyy = ER2.*sx./sy;
    ERzz = ER2.*sx.*sy;
    MURxx = ER2./sx.*sy;
    MURyy = ER2.*sx./sy;
    MURzz = ER2.*sx.*sy;

    if (verbose)
        toc();
        tic();
        println("(4) parsing from 2x to 1x grid");
    end

    # (4) parse down into 1x grid
    MURxx = MURxx[1:2:Nx2,2:2:Ny2];
    MURyy = MURyy[2:2:Nx2,1:2:Ny2];
    MURzz = MURzz[2:2:Nx2,2:2:Ny2];
    ERxx  = ERxx[2:2:Nx2,1:2:Ny2];
    ERyy  = ERyy[1:2:Nx2,2:2:Ny2];
    ERzz  = ERzz[1:2:Nx2,1:2:Ny2];
    (Nx,Ny) = size(ERxx);
    NGRID = [Nx,Ny];

    if (verbose)
        toc();
        tic();
        println("(5) constructing diagonal material matrices");
    end

    # (5) construct diagonal materials matrices
    ERxx_inv = spdiagm(1./ERxx[:]);
    #ERxx = spdiagm(ERxx[:]);
    ERyy_inv = spdiagm(1./ERyy[:]);
    #ERyy = spdiagm(ERyy[:]);
    #ERzz_inv = spdiagm(1./ERzz[:]);
    ERzz = spdiagm(ERzz[:]);
    MURxx_inv = spdiagm(1./MURxx[:]);
    #MURxx = spdiagm(MURxx[:]);
    MURyy_inv = spdiagm(1./MURyy[:]);
    #MURyy = spdiagm(MURyy[:]);
    #MURzz_inv = spdiagm(1./MURzz[:]);
    MURzz = spdiagm(MURzz[:]);

    if (verbose)
        toc();
        tic();
        println("(6) computing incident wave vector terms");
    end

    # (6) compute incident wave vector terms
    kinc = k0*n_ref*[sin(theta/180.0*pi); cos(theta/180.0*pi)];

    if (verbose)
        toc();
        tic();
        println("(7) calculating Derivative matrices");
    end

    # (7) calculate yee grid finite difference derivative matrices DEX,DEY,DHX,DHY
    (DEX,DEY,DHX,DHY) = yeeder(NGRID,k0*RES,BC,kinc/k0);

    if (verbose)
        toc();
        tic();
        println("(8) constructing system matrix A");
    end

    # (8) construct system matrix 'A' for the mode of interest
    if (Pol == "Hz")
        A = DHX*MURyy_inv*DEX + DHY*MURxx_inv*DEY+ERzz;
    elseif (Pol == "Ez")
        A = DEX*ERyy_inv*DHX + DEY*ERxx_inv*DHY+MURzz;
    else
        error = "throw error";
    end

    if (verbose)
        toc();
        tic();
        println("(9) computing source field without TFSF");
    end

    # (9) compute source field without TFSF
    f_src = zeros(Complex64,Nx,Ny);
    for i = (1:Nx)
        x = dx*i;
        for j = (1:Ny)
            y = dy*j;
            f_src[i,j] = exp(im*(kinc[1]*x + kinc[2]*y));
        end
    end

    if (verbose)
        toc();
        tic();
        println("(10) computing source");
    end

    # (10) calculate source
    if (true)
        # compute TFSF source vector b = (QA-AQ)f_src
        Q = spdiagm(Q[:]);
        b = (Q*A-A*Q)*f_src[:];
    end

    if (verbose)
        toc();
        tic();
        println("(11) solving Af=b for fields");
    end

    # (11) compute fields by solving f = inv(A)*b
    f = A\b;

    if (verbose)
        toc();
        tic();
        println("(12) converting into physical fields");
    end

    # (12) convert back into fields
    if (Pol == "Hz")
        Hz = reshape(f,Nx,Ny);
        Hx = Hz;
        Hy = Hz;
        Ex = Hz;
        Ey = Hz;
        Ez = Hz;
    else
        Ez = reshape(f,Nx,Ny);
        Ex = Ez;
        Ey = Ez;
        Hx = Ez;
        Hy = Ez;
        Hz = Ez;
    end

    if (verbose)
        toc();
        println("(~) finished with FDFD");
        toc();
    end

    return (Ex,Ey,Ez,Hx,Hy,Hz);
end

#BC = [0,0]
#Nx2 = 100;          Ny2 = 100;
#Nx = div(Nx2,2);     Ny = div(Ny2,2);
#ER2  = ones(Complex64,Nx2,Ny2);
#MUR2 = ones(Complex64,Nx2,Ny2);
#RES  = [1e-8,1e-8];
#NPML = [10 10 10 10]
#lambda0 = 1e-6;
#Pol = "Hz";
#theta = 10;
#Q = zeros(Int,Nx,Ny);
#Q[:,1:20] = 1;
#(Ex,Ey,Ez,Hx,Hy,Hz) = fdfd(ER2,MUR2,RES,NPML,BC,lambda0,Pol,theta,Q);
