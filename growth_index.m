function z = growth_index(len,lg)
    nel = length(len);
    z = 1;
    l = 0;
    for el = 1:nel
        l = l + len(nel+1-el);
        if nel-el>=1
            l1 = l + len(nel-el);
        end
        if l>=lg
            z = nel+1-el;
            break;
        elseif nel-el>=1 && l1>=lg && norm((l1 - lg))>=norm((l-lg))
                z = nel+1-el;
                break;
        end
    end
end