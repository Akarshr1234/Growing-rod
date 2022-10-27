function len = len_finder(pos,nel)
    len = zeros(nel,1);
    nnodes = length(pos)/3;
    npe = ((nnodes-1)/nel) + 1;
    ndof = npe*nnodes;
    L = zeros(npe*3,ndof,nel);

    for el = 1:nel
        L(:,(el-1)*(3*(npe-1)) + 1 :((el-1)*(3*(npe-1))) + npe*3,el) = eye(npe*3);
    end
    
    for el = 1:nel
        pos_el = L(:,:,el)*pos;
        pos_el1 = pos_el(1:3);
        pos_el2 = pos_el(end-2:end);
        len(el) = norm(pos_el1-pos_el2);
    end
end