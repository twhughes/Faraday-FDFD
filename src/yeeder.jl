function yeeder(NGRID,RES,BC,kinc=[0,0])

        #######################################################################
        # Calculate yee grid derivative operators
        #
        # Note: for normalized grid, use this function as follows:
        # (DEX,DEY,DHX,DHY) = yeeder(NGRID,k0*RES,BC,kinc/k0);
        #
        # NGRID = [Nx,Ny]   grid dimensions
        # RES = [dx,dy]     grid resolution of 1x grid
        # BC = [xbc,ybc]    boundary conditions
        #       bc = 0     -> Dirichlet (noes not require kinc)
        #       bc = -2    -> Periodic  (requires kinc)
        # kinc = [kx, ky]   incident wave vector (only required for periodic BC
        #
        #######################################################################


        (Nx,Ny) = NGRID;
        (dx,dy) = RES;
        (xbc,ybc) = BC;
        (kx,ky) = kinc;

        N = Nx*Ny;

        Lambda_x = dx*Nx;
        Lambda_y = dy*Ny;

        DEX = -1*speye(Complex64,N,N);
        DEY = -1*speye(Complex64,N,N);

        if (BC[1] == -1)
           kinc[1] = 0;
           BC[1] = -2;
        end
        if (BC[2] == -1)
           kinc[2] = 0;
           BC[2] = -2;
        end

        # Construct DE operators

        if (Nx > 1)
            DEX = DEX+spdiagm(ones(1,N)[:],1)[:,1:N];
        else
            DEX = DEX+spdiagm(zeros(1,N)[:],1)[:,1:N];
        end

        if (Ny > 1)
            DEY = DEY+spdiagm(ones(1,N)[:],Nx)[:,1:N];
        else
            DEY = DEY+spdiagm(zeros(1,N)[:],Nx)[:,1:N];
        end

        # fix DEX boundary conditions for Dirichlet boundary conditions
        for i in (Nx:Nx:N-Nx)
            DEX[i,i+1] = 0;
        end

        # if periodic boundary conditions in x

        if (xbc == -2)
            for i in (Nx:Nx:N)
               #Bloch DEX[i,i-Nx+1] = exp(im*Lambda_x*kx);
               DEX[i,i-Nx+1] = 1;#exp(im*Lambda_x*kx);
            end
        end

        if (ybc == -2)
           # for i in (N-Nx:N)
           #    DEY[i,i-N+Nx+1] = exp(im*Lambda_y*ky);
           # end
    #Bloch      #  DEY = DEY+transpose(spdiagm(exp(im*Lambda_y*ky)*ones(1,N)[:], N-Nx)[:,1:N]);
            DEY = DEY+transpose(spdiagm(1*ones(1,N)[:], N-Nx)[:,1:N]);
        end

        if (Nx == 1 && xbc == -2)
            DEX = (im*kx).*speye(N);
        end
        if (Ny == 1 && ybc == -2)
            DEY = (im*ky).*speye(N);
        end

        #added this 2 to fix problems
        DEX = DEX./dx;
        DEY = DEY./dy;

        DHX = -DEX';
        DHY = -DEY';
        return (DEX,DEY,DHX,DHY);
end

NGRID = [3,3];
BC = [0 0];
RES = [0.2,0.3];

(DEX,DEY,DHX,DHY) = yeeder(NGRID,RES,BC);
full(DEX);

