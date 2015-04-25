global q_el Na Ang joulesPerCalorie E_0 E_inf kB conv_factor;
q_el = 1.60217646e-19;
Na = 6.0221415e23;
Ang = 1e-10;
joulesPerCalorie = 4.184;
E_0 = 8.854187818e-12;
kB = Na * 1.3806488e-23/joulesPerCalorie/1000;
% now kB is in kcal/mol/K
T  = 300.0; % temp

strict_conv_factor = (Na/1000) * (q_el^2 / E_0) * 1e10 / 4 / pi / joulesPerCalorie;

conv_factor = 332.112; % this is what FFTSVD and other publications
                       % use, so we use it here.