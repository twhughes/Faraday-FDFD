function calc_gradient(Ex,lam,dl,beta)
    (Nx,Ny) = size(Ex);
 #   Np = 10;
 #   phis = (1:Np)*2*pi/Np;
 #   Gs = zeros(Np,1);
 #   G_max = 0.0;
 #   i_max = 0;
 #   E_tilde = zeros(Nx,1);
    xs = dl*(1:Nx)
 #   for pi in (1:Np)
 #       phi_1 = phis[pi];
 #       Gi = 0.0;
 #       for xi in (1:Nx)
 #           x = (xi-1)*dl;
 #           E_tilde[xi,1] = real( Ex[xi,Int(Ny/2) ]*exp(im*(2*pi/lam/beta*x+phi_1) ));
 #           Gi = Gi+real( Ex[xi,Int(Ny/2) ]*exp(im*(2*pi/lam/beta*x+phi_1) ));
 #       end
 #       Gs[pi] = Gi;
 #       if (Gi > G_max)
 #           G_max = Gi;
 #           i_max = pi;
 #       end
 #    end

    R = sum(real(Ex0[:,Int(Ny/2)].*exp(im*(2*pi/lam/beta*xs + 1*ones(Nx,1)))));
    I = sum(imag(Ex0[:,Int(Ny/2)].*exp(im*(2*pi/lam/beta*xs + 1*ones(Nx,1)))));

    G_max = abs(R+I)/Nx;
    phi_max = atan(abs(I)/abs(R));

    return (G_max,phi_max);
end
