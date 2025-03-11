data:=proc()
    p:=2^31-1:
    num_var:=3;
    str_var:=[seq(cat("x_",i),i=1..num_var)];# generates variable "x_{i}" as strings.
    vars:={seq(convert(str_var[i],name),i=1..num_var)};
    num_deg:=3:
    denom_deg:=3:
    # f:=randpoly(vars,sparse,degree=num_deg):
    # g:=randpoly(vars,sparse,degree=denom_deg):
    # B:=Construct_Blackbox(f,g,vars):
    # return B,vars,num_var,p:
end proc: