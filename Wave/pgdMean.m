function [amean] = pgdMean(pgd,a,threshold)
    ind = find(pgd>=threshold);
    amean = mean(a(ind),"all");
end

