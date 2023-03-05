function shape = SigmaSC(x, b)

phi_vec = 0:0.01:2*pi;
shape = zeros(2, length(phi_vec));
nr_Fourier_coeff = size(b,1);
for i = 1 : length(phi_vec)
    phi = phi_vec(i);
    R = calcFourierCoeff(phi, nr_Fourier_coeff);
    e = [cos(phi) sin(phi)]';
    shape(:, i) = R * b * e + x;
end

function fourie_coff = calcFourierCoeff(theta, nr_Fourier_coeff)

fourie_coff(1) = 1/2;

index = 1;
for i = 2 : 2 : (nr_Fourier_coeff - 1)
    fourie_coff(i : i + 1) = [cos(index * theta) sin(index * theta)];
    index = index + 1;
end