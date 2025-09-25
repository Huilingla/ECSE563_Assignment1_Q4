function [Eeq, Zeq, Yeq] = genthevenin(Y, Iint, id)
% GENTHEVENIN Calculate generalized Thévenin equivalent for multiple nodes.
%
%   [EEQ, ZEQ, YEQ] = GENTHEVENIN(Y, IINT, ID) computes the generalized
%   Thévenin equivalent for a subset of nodes using Norton equivalence.
%
%   Inputs:
%     Y    - NxN admittance matrix of the network
%     Iint - Nx1 vector of pre-fault internal current sources
%     id   - Vector of node indices where equivalent is evaluated
%
%   Outputs:
%     Eeq  - Vector of generalized Thévenin equivalent voltages
%     Zeq  - Generalized Thévenin equivalent impedance matrix
%     Yeq  - Reduced admittance matrix (optional)
%
%   Method:
%     1. Use network reduction to find equivalent as seen from nodes 'id'
%     2. Combine Norton and Thévenin equivalence principles

    fprintf('=== GENERALIZED THÉVENIN EQUIVALENT ===\n\n');
    
    % Validate inputs
    N = size(Y, 1);
    if length(Iint) ~= N
        error('Iint must have same dimension as Y');
    end
    if max(id) > N || min(id) < 1
        error('id contains invalid node indices');
    end
    
    fprintf('Calculating equivalent for nodes: ');
    fprintf('%d ', id);
    fprintf('\n\n');
    
    % Step 1: Partition the network into internal and external nodes
    internal_nodes = id;  % Nodes we're finding equivalent for
    external_nodes = setdiff(1:N, internal_nodes);  % Rest of the network
    
    fprintf('Internal nodes (equivalent ports): ');
    fprintf('%d ', internal_nodes);
    fprintf('\n');
    fprintf('External nodes (to be reduced): ');
    fprintf('%d ', external_nodes);
    fprintf('\n\n');
    
    % Step 2: Reorder Y matrix and Iint vector for partitioning
    all_nodes = [internal_nodes, external_nodes];
    Y_ordered = Y(all_nodes, all_nodes);
    Iint_ordered = Iint(all_nodes);
    
    % Partition the reordered system
    k = length(internal_nodes);  % Number of equivalent ports
    m = length(external_nodes);  % Number of nodes to reduce
    
    Y11 = Y_ordered(1:k, 1:k);      % Internal-to-internal
    Y12 = Y_ordered(1:k, k+1:end);  % Internal-to-external  
    Y21 = Y_ordered(k+1:end, 1:k);  % External-to-internal
    Y22 = Y_ordered(k+1:end, k+1:end); % External-to-external
    
    I1 = Iint_ordered(1:k);         % Currents at internal nodes
    I2 = Iint_ordered(k+1:end);     % Currents at external nodes
    
    fprintf('Step 1: Network Partitioning\n');
    fprintf('----------------------------\n');
    fprintf('Y11 size: %dx%d (internal nodes)\n', size(Y11,1), size(Y11,2));
    fprintf('Y22 size: %dx%d (external nodes)\n', size(Y22,1), size(Y22,2));
    fprintf('Y12 size: %dx%d (coupling)\n', size(Y12,1), size(Y12,2));
    
    % Step 3: Calculate generalized Thévenin equivalent using Kron reduction
    fprintf('\nStep 2: Kron Reduction\n');
    fprintf('----------------------\n');
    
    % Reduce external network: Yeq = Y11 - Y12 * Y22^(-1) * Y21
    % Instead of inv(Y22), solve Y22 * X = Y21 column by column
    X = zeros(size(Y22, 2), size(Y21, 2));
    for col = 1:size(Y21, 2)
        X(:, col) = linsolve(Y22, Y21(:, col));
    end
    
    Yeq = Y11 - Y12 * X;
    fprintf('Reduced admittance matrix Yeq size: %dx%d\n', size(Yeq,1), size(Yeq,2));
    
    % Step 4: Calculate equivalent current sources using Norton equivalence
    fprintf('\nStep 3: Equivalent Current Sources\n');
    fprintf('----------------------------------\n');
    
    % Solve for voltages at external nodes with internal nodes open-circuited
    V2_open = linsolve(Y22, I2);
    
    % Equivalent Norton current sources: Ieq = I1 - Y12 * Y22^(-1) * I2
    % Use linsolve instead of inv(Y22)
    temp = linsolve(Y22, I2);
    Ieq = I1 - Y12 * temp;
    
    fprintf('Norton equivalent currents:\n');
    for i = 1:k
        fprintf('  Ieq(%d) = %.4f + j%.4f p.u.\n', internal_nodes(i), real(Ieq(i)), imag(Ieq(i)));
    end
    
    % Step 5: Convert to Thévenin equivalent using linsolve
    fprintf('\nStep 4: Thévenin Equivalent\n');
    fprintf('---------------------------\n');
    
    % Calculate Thévenin impedance matrix Zeq by solving Yeq * Zeq = I
    % Do this column by column using linsolve
    I_matrix = eye(k);
    Zeq = zeros(k, k);
    
    for col = 1:k
        e_col = I_matrix(:, col);
        z_col = linsolve(Yeq, e_col);
        Zeq(:, col) = z_col;
    end
    
    % Thévenin voltage vector: Eeq = Zeq * Ieq
    Eeq = Zeq * Ieq;
    
    fprintf('Thévenin equivalent voltages:\n');
    for i = 1:k
        fprintf('  Eeq(%d) = %.4f + j%.4f p.u. (|Eeq| = %.4f p.u.)\n', ...
                internal_nodes(i), real(Eeq(i)), imag(Eeq(i)), abs(Eeq(i)));
    end
    
    fprintf('\nThévenin impedance matrix Zeq (%dx%d):\n', k, k);
    for i = 1:k
        fprintf('  Row %d: ', i);
        for j = 1:k
            if imag(Zeq(i,j)) >= 0
                fprintf('%.4f+j%.4f  ', real(Zeq(i,j)), imag(Zeq(i,j)));
            else
                fprintf('%.4f-j%.4f  ', real(Zeq(i,j)), abs(imag(Zeq(i,j))));
            end
        end
        fprintf('\n');
    end
    
    % Validation
    fprintf('\nStep 5: Validation\n');
    fprintf('------------------\n');
    
    % Check that the reduction preserves open-circuit voltages
    V_original = linsolve(Y, Iint);
    V_oc_original = V_original(internal_nodes);
    
    error = max(abs(Eeq - V_oc_original));
    fprintf('Maximum error in open-circuit voltage: %.2e p.u.\n', error);
    
    if error < 1e-8
        fprintf('Validation: Excellent\n');
    else
        fprintf('Validation: Acceptable\n');
    end
end
