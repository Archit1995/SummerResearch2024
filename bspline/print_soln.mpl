print_soln:=proc(System,Vars,params)
	local unordered_soln,Soln,soln,var,i,v:
	Soln:=get_eqn(Sys,Vars):
	unordered_soln:=convert(Soln,list):
	soln:=reording(unordered_soln,nops(Sys)):
	lprint("soln ",soln);
end proc: