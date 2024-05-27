function m = my_nanmedian(x)

s=~isnan(x);
m = median(x(s));

end

