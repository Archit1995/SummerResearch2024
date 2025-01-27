f1:=randpoly(x,degree=50):
f1:=f1/lcoeff(f1);
r:=rand(1000..1000000);
points:=[seq(r(),1..degree(f1)+1)]:
y:=[seq(eval(f1,x=points[i]),i=1..degree(f1)+1)]:
finterp:=interp(points,y,x);
f1-finterp;