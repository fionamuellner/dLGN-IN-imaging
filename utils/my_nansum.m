function m = my_nansum(x)

s=~isnan(x);
m = sum(x(s));

end

