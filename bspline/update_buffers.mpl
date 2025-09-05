update_buffers:=proc(mqrfr_results)
    local F,G,lg,deg_num,deg_den,num_points_mqrfr,num_eval,den_eval:
    print("In update_buffers"):
    F:=[seq(mqrfr_results[i][1],i=1..nops(mqrfr_results))]:
    G:=[seq(mqrfr_results[i][2],i=1..nops(mqrfr_results))]:
    lg:=[seq(mqrfr_results[i][3],i=1..nops(mqrfr_results))]:
    print("F: ",F):
    print("G: ",G):
    # Collect the degrees of the numerators and denominators
    deg_num:=[seq(degree(F[i],x),i=1..nops(F))]:
    deg_den:=[seq(degree(G[i],x),i=1...nops(G))]:
    print("deg_num: ",deg_num):
    print("deg_den: ",deg_den):
    # Calculate the number of points needed for MQRFR
    # num_points_mqrfr:=[seq(deg_num[i]+deg_den[i]+2,i=1..nops(deg_num))]:
    num_points_mqrfr:=max(deg_num)+max(deg_den)+2:
    print("num_points_mqrfr: ",num_points_mqrfr):
    # Evaluate the numerators and denominators at x=1
    num_eval:=[[seq(eval(F[i],x=1),i=1..nops(F))]]:
    den_eval:=[[seq(eval(G[i],x=1),i=1..nops(G))]]:
    print("num_eval: ",num_eval):
    print("den_eval: ",den_eval):
    return F,G,lg,deg_num,deg_den,num_points_mqrfr,num_eval,den_eval:
    end proc: