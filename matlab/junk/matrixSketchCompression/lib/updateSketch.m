function [B] = updateSketch(l, Ai, B)
%function Bout = updateSketch(l, Ai, B)


m = length(Ai);


B_hasZeroRows = false;
indZeroRow = -1;
for j = 1 : l
    if 1 && all(B(j,:)==0)
        B_hasZeroRows = true;
        indZeroRow = j;
        break;
    end
end

if (B_hasZeroRows)
    B(indZeroRow, :) = Ai;
else
    [U, S, V] = svd(B);
    
    diagS = diag(S);
    
    indexS=ceil(find(diagS~=0,1,'last')/2);
    delta=diagS(indexS)^2;
    
    
    I_l = eye(l);
    if (l >= m)
        I_delta = I_l(:, 1 : m);
    else
        d = m - l;
        I_delta = [I_l, zeros(l, d)];
    end
    %I_delta = I_l .* delta;
    I_delta = I_delta .* delta;
    S_check = sqrt(max((S .^ 2) - I_delta, 0));
    %B = S_check * V;
    B = S_check * V';
    
    
    for j = 1 : l
        %if (sum(B(j, :)) == 0)
        if 1 && all(B(j,:)==0)
            B_hasZeroRows = true;
            indZeroRow = j;
            break;
        end
    end
    B(indZeroRow, :) = Ai;
   
end




