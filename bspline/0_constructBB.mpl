Constuct_Sys_Blackbox:=proc(System,Vars,params) 
	local Lin_BB:
	Lin_BB:=proc(point_,p)option remember:
		local unordered_soln,Soln,soln,var,i,v:
		print("point_",point_):
		var:=params:
		Soln:=get_eqn(Sys,Vars):
		# Soln:=solve(Sys,Vars):
		# print("Soln ",Soln);
		unordered_soln:=convert(Soln,list):
		soln:=reording(unordered_soln,nops(Sys)):
		# print("soln ",soln);
		return [seq(eval(op(2,soln[i]),{seq(var[v]=point_[v],v=1..numelems(point_))}),i=1..nops(Sys))] mod p:
		# return [seq(eval(op(2,soln[i]),{seq(var[v]=point_[v],v=1..numelems(point_))}),i=1..nops(Sys))];
		# return [seq(0_constructBB)]
	end proc:
	return Lin_BB:
end proc:
# Solves the system of equations. 
get_eqn:=proc(Sys,vars)option remember:
	print("in get_eqn"):
	return solve(Sys,Vars):
end proc:
# extract_eqn:=proc(LBB)
 
# end proc:
