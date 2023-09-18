function [S_out] = mbinarize(S_in,option)
    [N,M] = size(S_in);
    option = 3;
    if sum(sum(isnan(S_in))) ~= 0
       S_in = ones(size(S_in)); 
    end
    if option == 1   
        th =(max(S_in(:))-min(S_in(:)))/2;
        S_in(S_in>=th)=1;%valores mayores que el umbral a 1
        S_in(S_in < th)=0; %menores que el umbral a 0
        S_out = S_in;
    elseif option == 2
        [idx,~] = kmeans(S_in(:),2);

        if idx(1)==2
            idx(idx==2) = 0;
        else
            idx(idx==1) = 0;
            idx(idx==2) = 1;
        end
        S_out = reshape(idx, [N M]);
    elseif option == 3
        S_out = S_in > 1e-5;
    end
end