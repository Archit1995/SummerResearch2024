with(LinearAlgebra);
 M := Matrix([
   [22, 44,  74,  1],
   [15, 14, -10, -2],
   [-25, -28, 20, 34]
 ]);
 R := ReducedRowEchelonForm(M);
R;
 M1 := copy(M):
 RowOperation(M1, 1, 1/22, inplace = true):
# R2 := R2 - 15*R1
RowOperation(M1, [2, 1], -15, inplace = true):

# R3 := R3 + 25*R1
RowOperation(M1, [3, 1], 25, inplace = true):

M1;
RowOperation(M1, 2, -1/16, inplace = true):

M1;
RowOperation(M1, [1, 2], -2, inplace = true):
RowOperation(M1, [3, 2], -22, inplace = true):

M1;
