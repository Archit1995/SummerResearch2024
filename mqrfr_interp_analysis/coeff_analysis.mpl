read "./2_mqrfr.mpl":
p:=19;
f:=3*x^2+4;
g:=x+2;
B:=f/g;

 g5:=5*x+2;
B5:=f/g5;

g3:=3*x+2;
B3:=f/g3;

xval:=[15,13,6,2,9];
m:=expand(product(x-xval[i],i=1..nops(xval)))mod p:

yval:=[seq(Eval(B,x=xval[i])mod p,i=1..nops(xval))];
u:=Interp(xval,yval,x)  mod p;

xval:=[15,13,6,2,1];
m:=expand(product(x-xval[i],i=1..nops(xval)))mod p:

yval:=[seq(Eval(B,x=xval[i])mod p,i=1..nops(xval))];
u:=Interp(xval,yval,x)  mod p;

xval:=[15,13,6,2,0];
m:=expand(product(x-xval[i],i=1..nops(xval)))mod p;

yval:=[seq(Eval(B,x=xval[i])mod p,i=1..nops(xval))]:
u:=Interp(xval,yval,x)  mod p;

xval:=[15,13,6,2,7];
m:=expand(product(x-xval[i],i=1..nops(xval)))mod p:

yval:=[seq(Eval(B,x=xval[i])mod p,i=1..nops(xval))];
u:=Interp(xval,yval,x)  mod p;

xval:=[15,13,6,2,5];
m:=expand(product(x-xval[i],i=1..nops(xval)))mod p:

yval:=[seq(Eval(B,x=xval[i])mod p,i=1..nops(xval))];
u:=Interp(xval,yval,x)  mod p;
# yval5:=[seq(Eval(B5,x=xval[i])mod p,i=1..nops(xval))];
# yval3:=[seq(Eval(B3,x=xval[i])mod p,i=1..nops(xval))];

# u:=Interp(xval,yval,x)  mod p;
# u5:=Interp(xval,yval5,x) mod p;
# u3:=Interp(xval,yval3,x) mod p;
# Quo(u,g,x) mod p;
# f,g,dq,lcg:=MQRFR(m,u,0,1,p);