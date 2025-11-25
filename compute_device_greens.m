function G_D_list = compute_device_greens(H_D, Sigma_L_list, Sigma_R_list, E_list, eta)
    N = length(E_list);
    dim = size(H_D,1);
    I = eye(dim);
    
    G_D_list = cell(N,1);

    fprintf('Computing device Green''s functions...\n');
    
    for idx = 1:N
        E = E_list(idx);
        Sigma_L = Sigma_L_list{idx};
        Sigma_R = Sigma_R_list{idx};

       
        if isempty(Sigma_L) || size(Sigma_L,1) ~= dim
            Sigma_L = zeros(dim,dim); 
        end
        if isempty(Sigma_R) || size(Sigma_R,1) ~= dim
            Sigma_R = zeros(dim,dim); 
        end

        A = (E + 1i*eta) * I - H_D - Sigma_L - Sigma_R;  
    
       
        try
            G_D = A \ I;
        catch
            fprintf('Warning: Backslash failed at E=%.3f, using pinv\n', E);
            G_D = pinv(A, 1e-10);
        end
        
        G_D_list{idx} = G_D;
        
        if mod(idx, 50) == 0
            fprintf('Completed %d/%d energies\n', idx, N);
        end
    end
end