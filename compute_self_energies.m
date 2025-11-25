% Authors:
% - Omid Soltani (omid.soltan@hotmail.com)
% - Mohammad Reza Jafari (mo.jafari@alzahra.ac.ir)
% - Aliasghar Shokri (aashokri@alzahra.ac.ir)
% 
% Alzahra University, Tehran, Iran
% Paper: Comparative transport characteristics of MXenes

function [Sigma_L_list, Sigma_R_list] = compute_self_energies(H0L, H0R, V, HLD, HDR, E_list, Vbias)
    eta = 1e-5;       
    max_iter = 200; 
    tol = 1e-8;

    E_F_shift = 0;  
    N = length(E_list);
    Sigma_L_list = cell(N,1);
    Sigma_R_list = cell(N,1);  

    for idx = 1:N
        E = E_list(idx);
        E_L = E + Vbias/2 + E_F_shift;  
        E_R = E - Vbias/2 + E_F_shift;  

        gL = surface_gf_sancho(H0L, V, E_L, eta, max_iter, tol);
        gR = surface_gf_sancho(H0R, V, E_R, eta, max_iter, tol);


     
        if isempty(gL) || isempty(gR) || any(isnan(gL(:))) || any(isnan(gR(:)))
            Sigma_L_list{idx} = zeros(size(HLD,2)); 
            Sigma_R_list{idx} = zeros(size(HDR,1)); 
            continue;
        end

   
        Sigma_L_list{idx} = HLD' * gL * HLD;
        Sigma_R_list{idx} = HDR * gR * HDR'; 
    end
end

