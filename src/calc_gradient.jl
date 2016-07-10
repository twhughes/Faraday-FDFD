function calc_gradient(Ex,lam,dl,beta)

    (Nx,Ny) = size(Ex);
    xs = dl*(1:Nx)
    Z = sum(Ex[:,Int(ceil(Ny/2))].*exp(im*(2*pi/lam/beta*xs)));
    R = real(Z);
    I = imag(Z);

    G_max = abs(Z)/Nx;
    phi_max = atan(abs(I)/abs(R));
    return (G_max,phi_max);
end
