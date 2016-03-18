function yeeder(NGRID,RES,BC,kinc)

        #######################################################################
        # Calculate yee grid derivative operators
        #
        # Note: for normalized grid, use this function as follows:
        # [DEX,DEY,DHX,DHY] = yeeder(NGRID,k0*RES,BC,kinc/k0);
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


end

yeeder(NGRID,RES,BC,kinc)
