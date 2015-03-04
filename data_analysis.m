clear all
close all

%% data files
filenames = {
'energy_matrix-rho0.30-T3.0.dat'
'energy_matrix-rho0.70-T0.3.dat'
'energy_matrix-rho0.70-T0.5.dat'
'energy_matrix-rho0.70-T0.7.dat'
'energy_matrix-rho0.70-T0.9.dat'
'energy_matrix-rho0.70-T1.0.dat'
'energy_matrix-rho0.70-T1.1.dat'
'energy_matrix-rho0.70-T1.5.dat'
'energy_matrix-rho0.80-T0.4.dat'
'energy_matrix-rho0.80-T0.6.dat'
'energy_matrix-rho0.80-T0.8.dat'
'energy_matrix-rho0.80-T1.0.dat'
'energy_matrix-rho0.80-T1.1.dat'
'energy_matrix-rho0.80-T2.0.dat'
'energy_matrix-rho0.88-T0.2.dat'
'energy_matrix-rho0.88-T0.3.dat'
'energy_matrix-rho0.88-T0.4.dat'
'energy_matrix-rho0.88-T0.5.dat'
'energy_matrix-rho0.88-T0.6.dat'
'energy_matrix-rho0.88-T0.7.dat'
'energy_matrix-rho0.88-T0.8.dat'
'energy_matrix-rho0.88-T0.9.dat'
'energy_matrix-rho0.88-T1.0.dat'
'energy_matrix-rho0.88-T1.1.dat'
'energy_matrix-rho1.20-T0.5.dat'
};

%% parameters
N = length(filenames);
N_part = 864;
N_equilibriation = 200;


datablock_size = 100;

Kb = 1;


for i = 1:N
    
    %% load data
    density(i) = wstr2num(filenames{i}(18:21));
    T_init(i)  = wstr2num(filenames{i}(24:26));
    M = load(filenames{i});
    E_kin(:,i) = M(:,3);
    E_pot(:,i) = M(:,4);
    time(:,i)  = M(:,2);
    step(:,i)  = M(:,1);
    virial(:,i)= M(:,5);
    
    %% parameters
    N_timesteps = length(step(:,i));
    N_measurement_length = N_timesteps - N_equilibriation;
    measurement_interval = N_equilibriation:1:N_timesteps;
    N_blocks = floor((N_measurement_length - N_equilibriation) / datablock_size);
    
%% measurements results
    kinenergy(i) = mean(E_kin(measurement_interval, i));
    potenergy(i) = mean(E_pot(measurement_interval, i));
    totenergy(i) = kinenergy(i) + potenergy(i);
    temperature(i) = 2*N_part/(3*(N_part-1)*Kb)*kinenergy(i);
    compressure(i) = 1 + 1/(3*Kb*temperature(i)*N_part)* mean(virial(measurement_interval, i));
    
%% time correlation length
    input = E_pot(measurement_interval,i);
    timecorrelation = xcorr(input - mean(input));
    zero_timecorrelation = length(input);
    tau_pot(i) = 1/2*sum(timecorrelation/timecorrelation(zero_timecorrelation));
    
    input = E_kin(measurement_interval,i);
    timecorrelation = xcorr(input - mean(input));
    zero_timecorrelation = length(input);
    tau_kin(i) = 1/2*sum(timecorrelation/timecorrelation(zero_timecorrelation));
    
%% data blocking
    for j = 1:N_blocks
        E_kin_block(j) = mean( E_kin((N_equilibriation+1+(j-1)*datablock_size):1:(N_equilibriation+1+j*datablock_size),i));
        E_pot_block(j) = mean( E_pot((N_equilibriation+1+(j-1)*datablock_size):1:(N_equilibriation+1+j*datablock_size),i));
        %E_kin_squared_block(j) = mean( (E_kin((N_equilibriation+1+(j-1)*datablock_size):1:(N_equilibriation+1+j*datablock_size),i)).^2);
    end 
    E_tot_block = E_pot_block + E_kin_block;
    error_E_kin(i) = std(E_kin_block);
    error_E_pot(i) = std(E_pot_block);
    error_E_tot(i) = std(E_tot_block);
    %error_E_kin_squared(i) = std(E_kin_squared_block)/sqrt(N_blocks);
    
    %% specific heat
    fluctuation_E_kin_squared(i) = mean((E_kin(:,i)-mean(E_kin(:,i))).^2);
    specificheat(i) = -1/(fluctuation_E_kin_squared(i)/(kinenergy(i))^2 - 2/(3*N_part));
    
end

%% plotting

figure
hold on
plot(error_E_kin)
plot(error_E_pot)
plot(error_E_tot)
title('error estimation energy')

figure % constant density
hold on
selected_indices = density == 0.88;
plot(temperature(selected_indices), compressure(selected_indices))
selected_indices = density == 0.80;
plot(temperature(selected_indices), compressure(selected_indices))
selected_indices = density == 0.70;
plot(temperature(selected_indices), compressure(selected_indices))
title('P(T)')
xlabel('Temperature')
ylabel('Compressibility factor')

figure % constant temperature
selected_indices = T_init == 1;
plot(density(selected_indices), compressure(selected_indices))
title('P(rho), T =~ 1.0')
xlabel('density')
ylabel('Compressibility factor')

figure %specific heat
selected_indices = density == 0.80;
plot(temperature(selected_indices), specificheat(selected_indices))
title('specificheat, rho =~ 0.88')
xlabel('Temperature')
ylabel('Compressibility factor')


%%

for i = 1:15
figure %energy fluctuations
hold on
plot(step(:,i), E_kin(:,i), '-r')
plot(step(:,i), E_pot(:,i) + E_kin(:,i), '-k')
plot(step(:,i), E_pot(:,i), '-b')

title(['Energy per particle, \rho = ', num2str(density(i)), ', T = ', num2str(T_init(i))], 'fontsize', 18)
xlabel('Timestep', 'fontsize', 16)
ylabel('Energy (\epsilon)', 'fontsize', 16)
legend('Kinetic energy', 'Total energy', 'Potential energy')
set(gca,'fontsize',14)
axis([0, 2501, -8, 4])
end


%%
figure
hold on 
plot(tau_kin)
plot(tau_pot)
title('correlation time')
