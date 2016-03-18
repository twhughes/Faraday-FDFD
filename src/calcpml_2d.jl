function calcpml_2d(NGRID,NPML)

        ###############################################################################################
        # Calculate UPML parameters sx,sy
        # NGRID = [Nx,Ny] grid dimensions
        # NPML = [Nxlo,Nxhi,Nylo,Nyhi] pml dimensions (low -> low index (1), high -> high index (Nx,Ny)
        ###############################################################################################

        eta0 = 376.73;     # free space impedence
        a_max = 3.0;       # 0 <= a_max <= 5
        p = 3.0;           # 3 <= p <= 5
        sigma_p_max = 1.0; # sigma_prime ~ 1

        (Nx,Ny) = NGRID;
        (Nxlo,Nxhi,Nylo,Nyhi) = NPML;

        sx = ones(Complex64,Nx,Ny);
        sy = ones(Complex64,Nx,Ny);

        for nx in (1:Nxlo)
                xp = nx/Nxlo;
                sig_p_x = sigma_p_max*(sin(pi*xp/2)^2);
                a = 1 + a_max*(xp)^p;
                sx[Nxlo-nx+1,:] = a*(1+im*eta0*sig_p_x);
        end
        for nx in (1:Nxhi)
                xp = nx/Nxhi;
                sig_p_x = sigma_p_max*(sin(pi*xp/2)^2);
                a = 1 + a_max*(xp)^p;
                sx[Nx-Nxhi+nx,:] = a*(1+im*eta0*sig_p_x);
        end
        for ny in (1:Nylo)
                yp = ny/Nylo;
                sig_p_y = sigma_p_max*(sin(pi*yp/2)^2);
                a = 1 + a_max*(yp)^p;
                sy[:,Nylo-ny+1] = a*(1+im*eta0*sig_p_y);
        end
        for ny in (1:Nyhi)
                yp = ny/Nyhi;
                sig_p_y = sigma_p_max*(sin(pi*yp/2)^2);
                a = 1 + a_max*(yp)^p;
                sy[:,Ny-Nyhi+ny] = a*(1+im*eta0*sig_p_y);
        end
        return (sx,sy);
end
