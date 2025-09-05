get_u:=proc(M,col,alpha,p)
	print("In get_u");
	F:=[seq(convert(M[..,i],list),i=1..col)];#Taking columns of M as list of lists.
	# print("F: ",F);
	U:=[seq(Interp(alpha,F[i],x) mod p,i=1..col)];
	# print("U: ",U);
	# print("__________________________________");
	return U;		
end proc;

