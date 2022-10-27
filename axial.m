function ax = axial(S)
%     axial gives the axial vector corresponding to a skew tensor S.
    S = skew(S);
    ax = zeros(3,1);
    ax(1) = S(3,2);
    ax(2) = S(1,3);
    ax(3) = S(2,1);
end
    