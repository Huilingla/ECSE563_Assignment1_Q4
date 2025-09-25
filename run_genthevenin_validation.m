% RUN_GENTHEVENIN_VALIDATION Validate generalized Thévenin equivalent
clear all; close all; clc;

% Load the IEEE 9-bus system data
ieee9_A1;

fprintf('=== GENERALIZED THÉVENIN EQUIVALENT VALIDATION ===\n\n');

% Calculate admittance matrix using linsolve-based approach
fprintf('Calculating admittance matrix...\n');
Y = admittance(nfrom, nto, r, x, b);

fprintf('IEEE 9-Bus System - Generalized Thévenin Analysis\n');
fprintf('=================================================\n\n');

%% Part (a): Nodes 3 and 5
fprintf('PART (A): NODES 3 AND 5\n');
fprintf('=======================\n\n');

id_a = [3, 5];
[Eeq_a, Zeq_a] = genthevenin(Y, Iint, id_a);

% Display results in a clear format
fprintf('\nSUMMARY - Nodes 3 and 5:\n');
fprintf('------------------------\n');
fprintf('Thévenin Voltages:\n');
fprintf('  Eeq(3) = %.4f ∠ %.2f° p.u.\n', abs(Eeq_a(1)), angle(Eeq_a(1))*180/pi);
fprintf('  Eeq(5) = %.4f ∠ %.2f° p.u.\n', abs(Eeq_a(2)), angle(Eeq_a(2))*180/pi);

fprintf('\nThévenin Impedance Matrix:\n');
fprintf('Zeq = [\n');
for i = 1:size(Zeq_a,1)
    fprintf('  ');
    for j = 1:size(Zeq_a,2)
        if imag(Zeq_a(i,j)) >= 0
            fprintf('%.4f+j%.4f  ', real(Zeq_a(i,j)), imag(Zeq_a(i,j)));
        else
            fprintf('%.4f-j%.4f  ', real(Zeq_a(i,j)), abs(imag(Zeq_a(i,j))));
        end
    end
    fprintf('\n');
end
fprintf('] p.u.\n');

%% Part (b): Nodes 9 and 4  
fprintf('\n\nPART (B): NODES 9 AND 4\n');
fprintf('=======================\n\n');

id_b = [9, 4];
[Eeq_b, Zeq_b] = genthevenin(Y, Iint, id_b);

% Display results
fprintf('\nSUMMARY - Nodes 9 and 4:\n');
fprintf('------------------------\n');
fprintf('Thévenin Voltages:\n');
fprintf('  Eeq(9) = %.4f ∠ %.2f° p.u.\n', abs(Eeq_b(1)), angle(Eeq_b(1))*180/pi);
fprintf('  Eeq(4) = %.4f ∠ %.2f° p.u.\n', abs(Eeq_b(2)), angle(Eeq_b(2))*180/pi);

fprintf('\nThévenin Impedance Matrix:\n');
fprintf('Zeq = [\n');
for i = 1:size(Zeq_b,1)
    fprintf('  ');
    for j = 1:size(Zeq_b,2)
        if imag(Zeq_b(i,j)) >= 0
            fprintf('%.4f+j%.4f  ', real(Zeq_b(i,j)), imag(Zeq_b(i,j)));
        else
            fprintf('%.4f-j%.4f  ', real(Zeq_b(i,j)), abs(imag(Zeq_b(i,j))));
        end
    end
    fprintf('\n');
end
fprintf('] p.u.\n');

%% Equivalent Circuit Sketch Information
fprintf('\n\nEQUIVALENT CIRCUIT FOR NODES 9 AND 4:\n');
fprintf('====================================\n\n');

fprintf('The equivalent circuit consists of:\n');
fprintf('1. Two voltage sources:\n');
fprintf('   - E9 = %.4f ∠ %.2f° p.u. at Node 9\n', abs(Eeq_b(1)), angle(Eeq_b(1))*180/pi);
fprintf('   - E4 = %.4f ∠ %.2f° p.u. at Node 4\n', abs(Eeq_b(2)), angle(Eeq_b(2))*180/pi);
fprintf('\n2. Impedance network with parameters:\n');
fprintf('   - Self-impedance at Node 9: Z99 = %.4f+j%.4f p.u.\n', real(Zeq_b(1,1)), imag(Zeq_b(1,1)));
fprintf('   - Self-impedance at Node 4: Z44 = %.4f+j%.4f p.u.\n', real(Zeq_b(2,2)), imag(Zeq_b(2,2)));
fprintf('   - Mutual impedance: Z94 = Z49 = %.4f+j%.4f p.u.\n', real(Zeq_b(1,2)), imag(Zeq_b(1,2)));

fprintf('\n3. Circuit topology:\n');
fprintf('   [E9]---[Z99]---[Z94]---[Z44]---[E4]\n');
fprintf('     |                 |\n');
fprintf('    Node 9           Node 4\n');

%% Validation with original system using linsolve
fprintf('\n\nVALIDATION WITH ORIGINAL SYSTEM:\n');
fprintf('==============================\n\n');

% Calculate original voltages using linsolve
V_original = linsolve(Y, Iint);

fprintf('Open-circuit voltages from original system:\n');
fprintf('Node  Original V(p.u.)  Thévenin Eeq(p.u.)  Error\n');
fprintf('----  ----------------  ------------------  ------\n');

for i = 1:length(id_b)
    node = id_b(i);
    orig_voltage = abs(V_original(node));
    thev_voltage = abs(Eeq_b(i));
    error = abs(orig_voltage - thev_voltage);
    fprintf('%2d      %8.4f           %8.4f        %.2e\n', ...
            node, orig_voltage, thev_voltage, error);
end

% Additional validation: Check Yeq * Zeq ≈ I
fprintf('\nMatrix Validation: Yeq * Zeq should ≈ Identity\n');
I_matrix = eye(size(Zeq_b,1));
Yeq_Zeq = Yeq * Zeq_b;
max_error = max(max(abs(Yeq_Zeq - I_matrix)));
fprintf('Maximum error in Yeq * Zeq = I: %.2e\n', max_error);

if max_error < 1e-10
    fprintf('Matrix inversion validation: Excellent\n');
else
    fprintf('Matrix inversion validation: Acceptable\n');
end

fprintf('\n=== GENERALIZED THÉVENIN ANALYSIS COMPLETE ===\n');
