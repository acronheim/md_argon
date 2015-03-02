N_equilibriation = 200;

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

N = length(filenames)
for i = 1:N
    density(i) = wstr2num(filenames{i}(18:21));
    T_init(i)  = wstr2num(filenames{i}(24:26));
    M = load(filenames{i});
    E_kin(:,i) = M(:,3);
    E_pot(:,i) = M(:,4);
    time(:,i)  = M(:,2);
    step(:,i)  = M(:,1);
    virial(:,i)= M(:,5);
end
