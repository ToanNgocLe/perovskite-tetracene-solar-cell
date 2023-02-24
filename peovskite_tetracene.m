clear;
clc;
resolution = 300;
h = 6.62607015 * 10^(-34);
c = 3* 10^8;
k_b = 1.380649 * 10^(-23);
Temp = 5672; 
G = zeros(1,resolution);
n = 1;
shortest_wavelength = 0.2;
longest_wavelength = 4;
increment = ((longest_wavelength-shortest_wavelength)/resolution);
lambda = zeros(1,resolution);
time_step = 200;

%parameters for perovskite
T_pev_base = 300; %K, base temperature of perovskite
T_pev_control = zeros(1,time_step);
T_pev_exp = zeros(1,time_step);
k_pev = 0.3; %W/(m*K)

%paremeters for tetracene
T_tet_base = 300;
T_tet = zeros(1,time_step);
k_tet = 0.16;


%bandgap of perovskite
bandgap_pev_base = 1.56; %eV
bandgap_pev_control = zeros(1,time_step);
bandgap_pev_control(1,1) = bandgap_pev_base;

bandgap_pev_exp = zeros(1,time_step);
bandgap_pev_exp(1,1) = bandgap_pev_base;


%thermal properties of perovskite
beta_pev = 311; %J/(kg*K)
%1J = 1.602*10^(19)eV
thickness_pev = 400*10^(-9); %m
l = 1; %m
volume_pev = thickness_pev * l * l;
bulk_density_pev = 4286.4; %kg/m^3
mass_pev = bulk_density_pev * volume_pev;
binding_energy_pev = 50*10^(-3); %eV
C_prime_pev = mass_pev * beta_pev;

%thermal properties of tetracene
beta_tet = 1585;
thickness_tet = 400*10^(-9);
volume_tet = thickness_tet * l * l;
density_tet = 1.2; %g/cm^3
mass_tet = density_tet * 10^(6) * 10^(-3) * volume_tet;
bandgap_tet = 2.43; 
C_prime_tet = mass_tet*beta_tet;

%electron-phonon interaction
A_eff = 8.09*10^(-3); % eV
hw_eff = 5.87 * 10^(-3); % eV

%thermal expansion influence
alpha_v = 1.57*10^(-4); %1/K
bulk_modulus = 18.8; %GPa
bandgap_change_pressure_pev = -50*10^(-3); %eV/GPa



%document of wavelength
for i = 1:1:resolution
    lambda(1,i) = shortest_wavelength + (i-1)*increment;
end


g_inv = zeros(1,resolution);
%graph for irradiance
figure
graph_space = linspace(shortest_wavelength,longest_wavelength,resolution);
for i = 1:1:resolution
G(1,n) = 0.473*((10293)*(lambda(1,i)^(-5)))* 1/(exp(2.5/lambda(1,i))-1);
g_inv(1,n) = G(1,n)*increment;
n = n+1;
end
plot(lambda,G);
g_total = sum(g_inv,'all');
% # of photons
phi = zeros(1,resolution);
for i = 1:1:resolution
    phi(1,i) = (G(1,i)* increment*10^(-6)*lambda(1,i))/(h*c);
end

%energy of photons
e_photons = zeros(1,resolution);
for i=1:1:resolution
    e_photons(1,i) = ((h*c)/(10^(-6)*lambda(1,i)))/(1.602*10^(-19));
end


%perovskite control case
e_w_control = zeros(resolution,time_step);
sum_e_w_control = zeros(1,time_step);
below_eg_control = zeros(resolution,time_step);
sum_below_eg_control = zeros(1,time_step);
for i = 1:1:time_step
    for j = 1:1:resolution
        if (e_photons(1,j) > (bandgap_pev_control(1,i)+binding_energy_pev))
            e_w_control(i,j) = phi(1,j)*e_photons(1,j) - phi(1,j)*(bandgap_pev_control(1,i)+binding_energy_pev);
            sum_e_w_control(1,i) = sum_e_w_control(1,i) + e_w_control(i,j);
        end
        if (e_photons(1,j) < bandgap_pev_control(1,i) + binding_energy_pev)
            below_eg_control(i,j) = phi(1,j)*e_photons(1,j);
            sum_below_eg_control(1,i) = sum_below_eg_control(1,i) + below_eg_control(i,j);
        end
    end

%temperature rise because of thermalization
    if i == 1
    T_pev_control(1,i) = T_pev_base;
    else
    T_pev_control(1,i) = T_pev_control(1,i-1) + (10^(-10))*(sum_e_w_control(1,i))/(1.602*10^(19))/C_prime_pev;
    end
    
    if T_pev_control(1,i) >= 65 + 273.15
        T_pev_control(1,i) = 65+ 273.15;
    end
electron_phonon_interaction = ((A_eff)/(4*T_pev_control(1,i)))*(hw_eff/(k_b*T_pev_control(1,i))) * (1/(sinh(hw_eff/(2*k_b*T_pev_control(1,i)))));
thermal_expansion_pev = -alpha_v * bulk_modulus * bandgap_change_pressure_pev;

if i == 1
    bandgap_pev_control(1,i+1) = bandgap_pev_control(1,i) + (thermal_expansion_pev + electron_phonon_interaction)*(T_pev_control(1,i) - T_pev_base);
else
    bandgap_pev_control(1,i+1) = bandgap_pev_control(1,i) + (thermal_expansion_pev + electron_phonon_interaction)*(T_pev_control(1,i) - T_pev_control(1,i-1));
end
end

e_total = 0;
fraction_sum_below_eg_control = zeros(1,time_step);
fraction_sum_e_w_control = zeros(1,time_step);
for i = 1:1:resolution
    e_total = e_total + phi(1,i)*e_photons(1,i);
end
for i = 1:1:time_step
fraction_sum_below_eg_control(1,i)= sum_below_eg_control(1,i)/e_total;
fraction_sum_e_w_control(1,i) = sum_e_w_control(1,i)/e_total;
end


%perovskite exp case
e_w_exp = zeros(resolution,time_step);
sum_e_w_exp = zeros(1,time_step);
below_eg_exp = zeros(resolution,time_step);
sum_below_eg_exp = zeros(1,time_step);
e_w_tet = zeros(resolution,time_step);
sum_e_w_tet = zeros(1,time_step);
below_eg_tet = zeros(resolution,time_step);
sum_below_eg_tet = zeros(1,time_step);

for i = 1:1:time_step
    for j = 1:1:resolution
        if (e_photons(1,j) > (bandgap_pev_exp(1,i)+binding_energy_pev) && e_photons(1,j) < bandgap_tet)
            e_w_exp(i,j) = phi(1,j)*e_photons(1,j) - phi(1,j)*(bandgap_pev_exp(1,i)+binding_energy_pev);
            sum_e_w_exp(1,i) = sum_e_w_exp(1,i) + e_w_exp(i,j);
            below_eg_tet(i,j) = phi(1,j)*e_photons(1,j);
            sum_below_eg_tet(1,i) = sum_below_eg_tet(1,i) + below_eg_tet(i,j);
        end
        if (e_photons(1,j) < bandgap_pev_exp(1,i) + binding_energy_pev)
            below_eg_exp(i,j) = phi(1,j)*e_photons(1,j);
            sum_below_eg_exp(1,i) = sum_below_eg_exp(1,i) + below_eg_exp(i,j);
            below_eg_tet(i,j) = phi(1,j)*e_photons(1,j);
            sum_below_eg_tet(1,i) = sum_below_eg_tet(1,i) + below_eg_tet(i,j);
        end
        if (e_photons(1,j) > bandgap_tet)
            e_w_tet(1,j) = phi(1,j)*e_photons(1,j) - phi(1,j)*bandgap_tet;
            sum_e_w_tet(1,i) = sum_e_w_tet(1,i) + e_w_tet(1,j);
        end
    end
 
    %temperature rise because of thermalization
    if i == 1
    T_pev_exp(1,i) = T_pev_base;
    else
    T_pev_exp(1,i) = T_pev_exp(1,i-1) + (10^(-10))*((sum_e_w_exp(1,i))/(1.602*10^(19)))/(C_prime_pev);
    end

    if i == 1
    T_tet(1,i) = T_tet_base;
    else
    T_tet(1,i) = T_tet(1,i-1) + (10^(-10))*((sum_e_w_tet(1,i))/(1.602*10^(19)))/(C_prime_tet);
    end


    % heat transfer
    T_pev_exp(1,i) = ((k_tet/thickness_tet)*T_tet(1,i) + (k_pev/thickness_pev)*T_pev_exp(1,i))/((k_tet/thickness_tet)+(k_pev/thickness_pev));
    T_tet(1,i) = T_pev_exp(1,i);
    
    if T_pev_exp(1,i) >= 65 + 273.15
        T_pev_exp(1,i) = 65 + 273.15;
        T_tet(1,i) = T_pev_exp(1,i);
    end


    %bandgap recalculation
electron_phonon_interaction = ((A_eff)/(4*T_pev_control(1,i)))*(hw_eff/(k_b*T_pev_control(1,i))) * (1/(sinh(hw_eff/(2*k_b*T_pev_control(1,i)))));
thermal_expansion_pev = -alpha_v * bulk_modulus * bandgap_change_pressure_pev;

if i == 1
    bandgap_pev_exp(1,i+1) = bandgap_pev_exp(1,i) + (thermal_expansion_pev + electron_phonon_interaction)*(T_pev_exp(1,i) - T_pev_base);
else
    bandgap_pev_exp(1,i+1) = bandgap_pev_exp(1,i) + (thermal_expansion_pev + electron_phonon_interaction)*(T_pev_exp(1,i) - T_pev_exp(1,i-1));
end


end

    %fraction for pev_exp
    e_total = 0;
fraction_sum_below_eg_exp = zeros(1,time_step);
fraction_sum_e_w_exp = zeros(1,time_step);
for i = 1:1:resolution
    e_total = e_total + phi(1,i)*e_photons(1,i);
end
for i = 1:1:time_step
fraction_sum_below_eg_exp(1,i)= sum_below_eg_exp(1,i)/e_total;
fraction_sum_e_w_exp(1,i) = sum_e_w_exp(1,i)/e_total;
end

    %fraction for tet
    e_total = 0;
fraction_sum_below_tet = zeros(1,time_step);
fraction_sum_e_w_tet = zeros(1,time_step);
for i = 1:1:resolution
    e_total = e_total + phi(1,i)*e_photons(1,i);
end
for i = 1:1:time_step
fraction_sum_below_tet(1,i)= sum_below_eg_tet(1,i)/e_total;
fraction_sum_e_w_tet(1,i) = sum_e_w_tet(1,i)/e_total;
end


figure
    plot(T_pev_control,'-b')
    hold on
    plot(T_pev_exp,'--r')
    hold off
    legend('control','experimental','location','east')
    xlabel('time (s)')
    ylabel('Temperature (K)')

x = [fraction_sum_e_w_control(1,1)*100; fraction_sum_below_eg_control(1,1)*100];
y = [fraction_sum_e_w_control(1,1)*100 fraction_sum_e_w_exp(1,1)*100; fraction_sum_below_eg_control(1,1)*100 fraction_sum_below_eg_exp(1,1)*100];
z = [fraction_sum_e_w_tet(1,1)*100; fraction_sum_below_tet(1,1)*100];

figure
bar(y)
hold off

figure
bar(z)



