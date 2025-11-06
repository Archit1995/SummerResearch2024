H:=[op(H),E];
B := Matrix(num_coeff, num_coeff, shape = identity);
B:=<B|Vector([seq(0,i=1..num_coeff)])>;
B:=<B;H>;
with(IntegerRelations):
LLL(B);
