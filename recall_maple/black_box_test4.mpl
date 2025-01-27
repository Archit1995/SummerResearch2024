B:=proc(a,b,point_)
    local c;
    c:= a/b;
    return seq(eval(c,x=i),i=point_[1..numelems(point_)]);
end proc;
f:=randpoly(x,degree=3):
g:=randpoly(x,degree=3):
r:=rand(1000..100000):
points:=[seq(r(),1..degree(f,x)+1)];
y:=B(f,g,points);