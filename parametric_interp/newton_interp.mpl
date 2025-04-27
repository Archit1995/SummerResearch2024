# newton interpolation input [x1,...xn], [y1,..,yn], output f(x)

x_val:=[11,22,3,14,7];
y_val:=[a,b,6,7,2];
degree_:=nops(x_val)-1;
order_:=1;
Slope:=[]:
coeff_:=[]:
# newton_interp:=proc(x_val,y_val,degree_,order_):    
    # if order_>degree_ then
    #     return coeff_;
    # end if
    # while(order_<=degree_)do 
        for i from 1 to nops(x_val) do 
            # print("i=",i);
            # Slope:=[op(Slope),(y_val[i+1]-y_val[i])/(x_val[i+order_]-x_val[i])];
            Slope:=[op(Slope),i+1];
        end do;
        # diff_;
        # for j in diff_ do 
        # coeff_:=[op(coeff_),diff_[1]];
        # order_:=order_+1;
        # y_val:=diff_;
        # diff_:=[];
        # end do;
        # newton_interp(x_val,diff_,degree_,order_);
