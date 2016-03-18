function yeeder(NGRID,RES,BC,kinc)

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
        #######################################################################


        (Nx,Ny) = NGRID;
        (dx,dy) = RES;
        (xbc,ybc) = BC;
        (kx,ky) = kinc;

        N = Nx*Ny;
        DEX = spzeros(N,N);
        DEY = spzeros(N,N);

   # Construct DE operators
        # diagonal terms
        for i in (1:N)
            DEX[i,i] = -1;
            DEY[i,i] = -1
        end

        # off diagonal terms for DEX
        for i in (1:N-1)
            DEX[i,i+1] = 1;
        end

        # off diagonal terms for DEY
        for i in (1:N-Nx)
            DEY[i,i+Nx] = 1;
        end

        # fix DEX boundary conditions for Dirichlet boundary conditions
        if (xbc == 0)
            for i in (Nx:Nx:N-Nx)
                DEX[i,i+1] = 0;
            end
        end



        DEX = DEX./dx;
        DEY = DEY./dy;

        DHX = -DEX';
        DHY = -DEY';
        return (DEX,DEY,DHX,DHY)
end

NGRID = [3 2];
RES = [0.1 1];
BC = [0 0];
kinc = [0 0];
(DEX,DEY,DHX,DHY) = yeeder(NGRID,RES,BC,kinc);
full(DEX)
#full(DEY);
