#This example demonstrates that we cannot interpolate rational functions by 
#interpolating the numerator and denominator separately.
# Thus provides motivation for MQRFR. 
B:=proc(a,b,point_)
    local c;
    c:= a/b;
    return [seq(eval(c,x=i),i=point_[1..numelems(point_)])];
end proc;
extract_numerator:=proc(r)
    return op(1,r);
end proc:
extract_denominator:=proc(r)
    return op(2,r);
end proc:
f:=randpoly(x,degree=3);
f:=f/lcoeff(f);
g:=randpoly(x,degree=3);
g:=g/lcoeff(g);
r:=rand(1000..100000):
points:=[seq(r(),1..degree(f,x)+1)];
y:=B(f,g,points);
T_:=interp(points,y,x):
Num_:=[seq(extract_numerator(i),i=y[1..numelems(y)])];
Denom_:=[seq(extract_denominator(i),i=y[1..numelems(y)])];
N_:=interp(points,Num_,x);
D_:=interp(points,Denom_,x);
N_:=N_/lcoeff(N_);
D_:=D_/lcoeff(D_);
# f-N_;
# g-D_;
yf:=[seq(eval(f,x=i),i=points[1..numelems(points)])];
yg:=[seq(eval(g,x=i),i=points[1..numelems(points)])];
interp(points,yf,x);
interp(points,yg,x);
