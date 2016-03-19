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

        DEX = spzeros(Complex64,N,N);
        DEY = spzeros(Complex64,N,N);

   # Construct DE operators
        # diagonal terms
        for i in (1:N)
            DEX[i,i] = -1;
            DEY[i,i] = -1;
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
        for i in (Nx:Nx:N-Nx)
            DEX[i,i+1] = 0;
        end

        # if periodic boundary conditions in x

        if (xbc == -2)
            for i in (Nx:Nx:N)
               DEX[i,i-Nx+1] = exp(im*Lambda_x*kx);
            end
        end

        if (ybc == -2)
            for i in (N-Nx:N)
               DEY[i,i-N+Nx+1] = exp(im*Lambda_y*ky);
            end
        end

        if (Nx == 1)
            DEX = (im*kx).*speye(N);
        end
        if (Ny == 1)
            DEY = (im*ky).*speye(N);
        end

        DEX = DEX./dx;
        DEY = DEY./dy;

        DHX = -DEX';
        DHY = -DEY';
        return (DEX,DEY,DHX,DHY);
end
