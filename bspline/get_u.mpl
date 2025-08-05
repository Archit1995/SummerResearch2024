get_u:=proc(M,col,alpha,p)
	F:=[seq(convert(M[..,i],list),i=1..col)];
	U:=[seq(Interp(alpha,F[i],x) mod p,i=1..col)];
	return U;		
end proc;