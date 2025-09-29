function Y = admittance(nfrom, nto, r, x, b)
% ADMITTANCE Calculate the admittance matrix for an AC power network using KCL.
%
%   Y = ADMITTANCE(NFROM, NTO, R, X, B) computes the nodal admittance matrix
%   by following the nodal analysis steps from the lecture:
%     1. Choose reference node
%     2. Use KCL at each node
%     3. Express branch currents in terms of node voltages: YV = I
%
%   Inputs:
%     nfrom - Mx1 vector of branch starting nodes
%     nto   - Mx1 vector of branch ending nodes  
%     r     - Mx1 vector of branch resistances (p.u.)
%     x     - Mx1 vector of branch reactances (p.u.)
%     b     - Mx1 vector of branch shunt susceptances (p.u.)
%
%   Output:
%     Y     - NxN nodal admittance matrix (p.u.)

    % Step 1: Determine network size and choose reference node
    all_nodes = unique([nfrom; nto]);
    N = max(all_nodes);  % Total number of nodes
    M = length(nfrom);   % Total number of branches
    
    % Validate input dimensions
    if length(nto) ~= M || length(r) ~= M || length(x) ~= M || length(b) ~= M
        error('All input vectors must have the same length');
    end
    
    fprintf('Step 1: Network has %d nodes, %d branches\n', N, M);
    fprintf('        Reference node chosen as node %d (ground)\n\n', N);
    
    % Step 2: Initialize admittance matrix
    Y = zeros(N, N);
    
    fprintf('Step 2: Building admittance matrix using KCL principles\n');
    fprintf('        Diagonal elements: sum of admittances connected to node\n');
    fprintf('        Off-diagonal elements: -admittance between nodes\n\n');
    
    % Step 3: Build admittance matrix by applying KCL at each node
    for m = 1:M
        i = nfrom(m);
        j = nto(m);
        
        % Calculate series admittance for this branch
        z_series = r(m) + 1i*x(m);
        y_series = 1/z_series;
        
        % Shunt admittance (half at each end) from line charging
        y_shunt = 1i * b(m) / 2;
        
        % Apply KCL principles:
        % - Current leaving node i due to branch i-j: y_series * (Vi - Vj)
        % - This contributes +y_series to Y(i,i) and -y_series to Y(i,j)
        
        % Update diagonal elements (sum of admittances connected to node)
        Y(i, i) = Y(i, i) + y_series + y_shunt;
        Y(j, j) = Y(j, j) + y_series + y_shunt;
        
        % Update off-diagonal elements (-admittance between nodes)
        Y(i, j) = Y(i, j) - y_series;
        Y(j, i) = Y(j, i) - y_series;
        
        fprintf('  Branch %d-%d: y_series = %.4f + j%.4f, y_shunt = j%.4f\n', ...
                i, j, real(y_series), imag(y_series), imag(y_shunt));
    end
    
    fprintf('\nStep 3: Admittance matrix Y constructed such that YV = I\n');
    fprintf('        where V is node voltage vector and I is current injection vector\n\n');
end
