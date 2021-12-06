function pos = generate_initial_population(D,ps,is_constrained_DOF,GSCP, basic_node, DOF_constrained, N_node,sym_M,benchmark)
disp('initialisation is in progress')
pos = randi([0,1],ps,D);
for k = 1: ps
    X = pos(k,:);
    offspring_pass = false;
    ctr = 0;
    while offspring_pass == false
       ctr = ctr + 1;
       kp = key_p(X);
       f_x = fitness_topology(X,is_constrained_DOF,GSCP, basic_node, DOF_constrained, N_node,sym_M,benchmark);
        if sum(f_x == 0) == 3 % feasible
            % visited before : no calculation
            offspring_pass = true;
            sprintf('%d individual is established',k)
        else % infeasible
            % repair
            [Y,ky] = ea_1p1_bitflip(X);
            f_y = fitness_topology(Y,is_constrained_DOF,GSCP, basic_node, DOF_constrained, N_node, sym_M, benchmark);
            if lex_compare(f_y,f_x)
                X = Y;
                f_x = f_y;
            end
        end
    end
    pos(k,:) = X;
end