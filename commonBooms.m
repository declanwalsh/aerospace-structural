% takes boom numbers from 2 cells
% finds common boom pairs (i.e. shared panels)
% only works for panels with no booms in between currently
function [ commonB, idx1, idx2 ] = commonBooms( B1, B2 )

    [commonB, i1, i2]  = intersect(B1, B2);
    
    if (length(i1) > 1)
        if (max(i1) == length(B1))
            idx1 = length(B1);
        else
            idx1 = min(i1);
        end
        if (max(i2) == length(B2))
            idx2 = length(B2);
        else
            idx2 = min(i2);
        end
    else
        fprintf('Only one common boom found - no panel')
        idx1 = 0;
        idx2 = 0;
    end
        
end

