construct_rational_function:=proc(f,g)
    return f/g;
end proc:
extract_numerator:=proc(r)
    return op(r,1);
end proc:
extract_denominator:=proc(r)
    return op(r,2);
end proc:
a:=randpoly(x,degree=3):
b:=randpoly(x,degree=5):
r:=construct_rational_function(a,b);
sigma1:=eval(r,x=1);
