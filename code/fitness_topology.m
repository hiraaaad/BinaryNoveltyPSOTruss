function fitness = fitness_topology(X,is_constrained_DOF,GSCP, basic_node, DOF_constrained, N_node,sym_M,benchmark)

switch benchmark
    case {10 15, 68}
        D = 2;
        [NAP,NAPconnector]=find_NAP(GSCP,X);
        grubler = D*sum(NAP)-sum(X)-numel(DOF_constrained);
    case {25, 72, 224}
        D = 3;
        M(sym_M(1,:)) = X;
        M(sym_M(2:end,:))= repmat(M(sym_M(1,:)),size(sym_M,1)-1,1);
        [NAP,NAPconnector]=find_NAP(GSCP,M);
        grubler = D*sum(NAP)-sum(M)-numel(DOF_constrained);
    case {52, 200, 47}
        D = 2;
        M(sym_M(1,:)) = X;
        M(sym_M(2:end,:))= repmat(M(sym_M(1,:)),size(sym_M,1)-1,1);
        [NAP,NAPconnector]=find_NAP(GSCP,M);
        grubler = D*sum(NAP)-sum(M)-numel(DOF_constrained);
end


fitness_basic_nodes = sum(NAP(basic_node)==0);

fitness_grubler = max(0,grubler); 

DOFonNAP = sum(reshape(is_constrained_DOF,D,N_node));
NAPconnector_all=NAPconnector+DOFonNAP-D;
fitness_equilibrium = abs(sum(min(NAPconnector_all(NAP==1),0)));

fitness = [fitness_basic_nodes,fitness_grubler,fitness_equilibrium];
end