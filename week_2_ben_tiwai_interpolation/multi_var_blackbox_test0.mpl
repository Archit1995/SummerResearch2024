with(LinearAlgebra):
with(ArrayTools):
# each point is a n dimensional vector, where n is the number of variables.
# B takes multiple n dimensional vectors as input and evaluates the polynomial at each vector.
# B outputs a scalar value corresponsing to each input vector.
B:=proc(a,var,point_)
    local u,v;
     # v:=[seq(eval(a,{var[v]=point_[v]}),v=1..numelems(var))]; This evaluates 1 var at t time
    # return eval(a,{seq(var[v]=point_[v],v=1..numelems(point_))});
    return [seq(eval(a,{seq(var[v]=point_[u][v],v=1..numelems(point_[u]))}),u=1..numelems(point_))];
end proc;
vars:={x,y,z}:
f:=randpoly(vars,degree=5);
deg:=degree(f):
points:=[[seq(i,i=1..3)],[seq(i,i=4..6)]];
B(f,vars,points);