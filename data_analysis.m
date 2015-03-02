%% data files
filenames = {
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
};

%% parameters
N = length(filenames);
N_part = 864;
N_equilibriation = 200;
N_timesteps = 2501;
N_measurement_length = N_timesteps - N_equilibriation;
measurement_interval = N_equilibriation:1:N_timesteps;
datablock_size = 100;
N_blocks = (N_measurement_length - N_equilibriation) / datablock_size;
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
    
%% measurements results
    kinenergy(i) = mean(E_kin(measurement_interval, i));
    potenergy(i) = mean(E_pot(measurement_interval, i));
    totenergy(i) = kinenergy(i) + potenergy(i);
    temperature(i) = 2*N_part/(3*(N_part-1)*Kb)*kinenergy(i);
    compressure(i) = 1 + 1/(3*Kb*Temperature(i)*N_part)* mean(virial(measurement_interval, i));
    
    %% data blocking
    for j = 1:N_blocks
        E_kin_block(j) = mean( E_kin((N_equilibriation+1+(j-1)*datablock_size):1:(N_equilibriation+1+j*datablock_size),i));
        E_pot_block(j) = mean( E_pot((N_equilibriation+1+(j-1)*datablock_size):1:(N_equilibriation+1+j*datablock_size),i));
        %E_kin_squared_block(j) = mean( (E_kin((N_equilibriation+1+(j-1)*datablock_size):1:(N_equilibriation+1+j*datablock_size),i)).^2);
    end 
    E_tot_block = E_pot_block + E_kin_block;
    error_E_kin(i) = std(E_kin_block)/sqrt(N_blocks);
    error_E_pot(i) = std(E_pot_block)/sqrt(N_blocks);
    error_E_tot(i) = std(E_tot_block)/sqrt(N_blocks);
    %error_E_kin_squared(i) = std(E_kin_squared_block)/sqrt(N_blocks);
    
    %% specific heat
    fluctuation_E_kin_squared(i) = mean((E_kin(:,i)-mean(E_kin(:,i))).^2);
    specificheat(i) = -1/(fluctuation_E_kin_squared(i)/(kinenergy(i))^2 - 2/(3*N_part));
    
end

    error_E_kin
    error_E_pot
    error_E_tot
%% plotting

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
figure %energy fluctuations
hold on
i = 16
plot(step(:,i), E_kin(:,i), '-r')
plot(step(:,i), E_pot(:,i), '-b')
plot(step(:,i), E_pot(:,i) + E_kin(:,i), '-k')
title('equilibriation')
xlabel('timestep')
ylabel('Energy')

