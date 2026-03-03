    |\^/|     Maple 2023 (X86 64 LINUX)
._|\|   |/|_. Copyright (c) Maplesoft, a division of Waterloo Maple Inc. 2023
 \  MAPLE  /  All rights reserved. Maple is a trademark of
 <____ ____>  Waterloo Maple Inc.
      |       Type ? for help.
> read "./0_construct_BB.mpl":
> read "./1_MRFI.mpl":
> read "./2a_NDSA.mpl":
> read "./2b_Deterministic_NDSA.mpl":
> read "./3_get_point_on_affine_line.mpl":
> read "./4_MQRFR.mpl":
> read "./5_BMEA.mpl":
> read "./6_generate_monomials.mpl":
> read "./7_zippel_vandermonde_solver.mpl":
> read "./8_construct_final_polynomial.mpl":
> read "./9_rvr.mpl":
> read "./data_gen.mpl":
> read "./helpers.mpl":



# Vars,F,G,num_vars,num_eqn,params := get_data(1):
# counter := 0:
# num_lines:=0:
# B:= Construct_Rational_Blackbox(F,G,Vars):



# test_case:="rand":
# num_var:=3:
# num_terms:=11:
# den_terms:=9:
# Vars,F,G,num_vars,num_eqn,params:=get_data(test_case,num_var,num_terms,den_terms):
# counter := 0:
# num_lines:=0:
# B:= Construct_Rational_Blackbox(F,G,Vars):


# test_case:="rat_rand":
# num_var:=3:
# num_terms:=11:
# den_terms:=9:
# num_coeff_bound:=20:
# den_coeff_bound:=20:
# Vars,F,G,num_vars,num_eqn,params:=get_data(test_case,num_var,num_terms,den_terms,num_coeff_bound,den_coeff_bound):
# counter := 0:
# B:= Construct_Rational_Blackbox(F,G,Vars):

# lprint("Variables:", Vars):
# lprint("Numerator F:", F):
# lprint("Denominator G:", G):
# print("num_eqn =",num_eqn):



# test_case:="example":
# test_case:="small_sys_low_deg":
#  test_case:="bspline":
# test_case:="small_Sys":
> test_case:="mike":
# test_case:="bsbug":
> num_lines:=0:
> Sys, Vars, params, num_vars, num_eqn:= get_data(test_case):
                                 "in get_data"

> counter := 0:
> B := Constuct_Sys_Blackbox(Sys, Vars, params):


# p:= 2^31 - 1:
> p:=107:



# print("Number of equations:", num_eqn):
# print("Number of parameters:", num_vars):
# Create black box
> try
> Num,Den := MRFI(B, num_vars, num_eqn, params, p):
>     catch:
>     lprint("ERROR:", lasterror()):
> end try:
"MRFI  ========================================"
"MRFI Starting MRFI"
"MRFI Number of parameters:", 2
"MRFI Number of equations:", 2
"MRFI ========================================"
"MRFI numerator_done: ", [false, false]
"MRFI denominator_done: ", [false, false]
                                   "In NDSA"

                                   "T:= ", 4

"NDSA: alpha: ", [2, 3, 4, 5]
"NDSA:m: ", x^4+93*x^3+71*x^2+60*x+13
"NDSA:Psi_alpha: ", [[2, 8], [3, 15], [4, 22], [5, 29]]
             "NDSA: Y = ", [[41, 13], [45, 98], [46, 88], [56, 94]]

"row: ", 4, " col: ", 2
          "---------------------------------------------------------"

                                   "In MQRFR"

"q[", 1, "]=", 54*x+38, " degree q=", 1
"f=", 2*x^3+34*x^2+10*x+83
"g=", 1
                                       3       2
          "NDSA: result", 1, ": ", [2 x  + 34 x  + 10 x + 83, 1, 1, 1]

          "---------------------------------------------------------"

                                   "In MQRFR"

"q[", 1, "]=", 55*x+86, " degree q=", 1
"f=", 72*x^3+x+77
"g=", 1
                                             3
               "NDSA: result", 2, ": ", [72 x  + x + 77, 1, 1, 1]

               "NDSA:MQRFR failed. Trying again with more points"

                     "NDSA: mqrfr_status: ", [false, false]

"NDSA: Resetting result for equation ", 1
"NDSA: Resetting result for equation ", 2
               "_______________________________________________"

                                   "T:= ", 8

"NDSA: alpha: ", [2, 3, 4, 5, 6, 7, 8, 9]
"NDSA:m: ", x^8+63*x^7+77*x^6+43*x^5+79*x^4+73*x^3+5*x^2+25*x+43
"NDSA:Psi_alpha: ", [[2, 8], [3, 15], [4, 22], [5, 29], [6, 36], [7, 43], [8, 
50], [9, 57]]
"NDSA: Y = ", [[41, 13], [45, 98], [46, 88], [56, 94], [42, 83], [82, 71],

    [49, 18], [86, 33]]

"row: ", 8, " col: ", 2
          "---------------------------------------------------------"

                                   "In MQRFR"

"q[", 1, "]=", 81*x+9, " degree q=", 1
"f=", 37*x^7+68*x^6+94*x^5+37*x^4+55*x^3+28*x^2+55*x+66
"g=", 1
"q[", 4, "]=", 40*x^4+38*x^3+46*x^2+84*x+103, " degree q=", 4
"f=", 49*x+49
"g=", 22*x^3+27*x^2+49*x
                                               3       2
         "NDSA: result", 1, ": ", [46 x + 46, x  + 45 x  + 46 x, 4, 22]

          "---------------------------------------------------------"

                                   "In MQRFR"

"q[", 1, "]=", 28*x+82, " degree q=", 1
"f=", 65*x^7+45*x^6+56*x^5+69*x^4+85*x^3+35*x^2+27*x+23
"g=", 1
"q[", 3, "]=", 99*x^5+62*x^4+88*x^3+62*x^2+46*x+94, " degree q=", 5
"f=", 19*x+88
"g=", 19*x^2+106*x+18
                                               2
           "NDSA: result", 2, ": ", [x + 106, x  + 45 x + 46, 5, 19]

                   "NDSA: MQRFR successful for equation ", 1

                   "NDSA: MQRFR successful for equation ", 2

                       "NDSA: Termination condition met"

                       "MRFI num_points_mqrfr: ", [6, 5]

                        "MRFI max_num_points_mqrfr: ", 6

  "______________________________________________________________________________"

                       "in main evaluation loop of MRFI"

"MRFI T_old=", 1
"MRFI T=", 4
                       "MRFI sigma_[", 1, "] = ", [2, 3]

                            "In Deterministic_NDSA"

               "Deterministic_NDSA: alpha: ", [2, 3, 4, 5, 6, 7]

"Deterministic_NDSA:Psi_alpha: ",

    [[2, 3], [3, 10], [4, 17], [5, 24], [6, 31], [7, 38]]

"Deterministic_NDSA: Y: ",

    [[69, 92], [15, 21], [95, 39], [49, 101], [31, 94], [13, 33]]

                              6       5       4       3       2
   "Deterministic_NDSA:m: ", x  + 80 x  + 81 x  + 47 x  + 75 x  + 104 x + 11

"Deterministic_NDSA: u: ", 20*x^5+72*x^4+51*x^3+6*x^2+105*x+96
          "---------------------------------------------------------"

                                   "In MQRFR"

"q[", 1, "]=", 91*x+83, " degree q=", 1
"f=", 20*x^5+72*x^4+51*x^3+6*x^2+105*x+96
"g=", 1
"q[", 4, "]=", 86*x^2+20*x+23, " degree q=", 2
"f=", 62*x+62
"g=", 6*x^3+67*x^2+62*x
                                   5       4       3       2
        "Deterministic_NDSA:m: ", x  + 87 x  + 48 x  + 62 x  + 81 x + 29

"Deterministic_NDSA: u: ", 15*x^4+50*x^3+37*x+20
          "---------------------------------------------------------"

                                   "In MQRFR"

"q[", 1, "]=", 50*x+46, " degree q=", 1
"f=", 15*x^4+50*x^3+37*x+20
"g=", 1
"q[", 3, "]=", 25*x^2+56*x+40, " degree q=", 2
"f=", 64*x+89
"g=", 64*x^2+37*x+55
  "______________________________________________________________________________"

                       "MRFI sigma_[", 2, "] = ", [4, 9]

                            "In Deterministic_NDSA"

               "Deterministic_NDSA: alpha: ", [2, 3, 4, 5, 6, 7]

"Deterministic_NDSA:Psi_alpha: ",

    [[2, 102], [3, 2], [4, 9], [5, 16], [6, 23], [7, 30]]

"Deterministic_NDSA: Y: ",

    [[89, 72], [97, 46], [21, 6], [23, 20], [24, 101], [20, 26]]

                              6       5       4       3       2
   "Deterministic_NDSA:m: ", x  + 80 x  + 81 x  + 47 x  + 75 x  + 104 x + 11

"Deterministic_NDSA: u: ", 49*x^5+22*x^4+14*x^3+65*x^2+60*x+31
          "---------------------------------------------------------"

                                   "In MQRFR"

"q[", 1, "]=", 83*x+67, " degree q=", 1
"f=", 49*x^5+22*x^4+14*x^3+65*x^2+60*x+31
"g=", 1
"q[", 4, "]=", 56*x^2+82*x+82, " degree q=", 2
"f=", 18*x+18
"g=", 19*x^3+86*x^2+18*x
                                   5       4       3       2
        "Deterministic_NDSA:m: ", x  + 87 x  + 48 x  + 62 x  + 81 x + 29

"Deterministic_NDSA: u: ", 20*x^4+88*x^3+27*x^2+77*x+70
          "---------------------------------------------------------"

                                   "In MQRFR"

"q[", 1, "]=", 91*x+48, " degree q=", 1
"f=", 20*x^4+88*x^3+27*x^2+77*x+70
"g=", 1
"q[", 3, "]=", x^2+62*x+82, " degree q=", 2
"f=", x+43
"g=", x^2+89*x+46
  "______________________________________________________________________________"

                       "MRFI sigma_[", 3, "] = ", [8, 27]

                            "In Deterministic_NDSA"

               "Deterministic_NDSA: alpha: ", [2, 3, 4, 5, 6, 7]

"Deterministic_NDSA:Psi_alpha: ",

    [[2, 92], [3, 99], [4, 106], [5, 6], [6, 13], [7, 20]]

"Deterministic_NDSA: Y: ",

    [[35, 19], [17, 19], [62, 72], [67, 83], [49, 76], [36, 10]]

                              6       5       4       3       2
   "Deterministic_NDSA:m: ", x  + 80 x  + 81 x  + 47 x  + 75 x  + 104 x + 11

"Deterministic_NDSA: u: ", 41*x^5+41*x^4+44*x^3+87*x^2+21*x
          "---------------------------------------------------------"

                                   "In MQRFR"

"q[", 1, "]=", 47*x+75, " degree q=", 1
"f=", 41*x^5+41*x^4+44*x^3+87*x^2+21*x
"g=", 1
"q[", 4, "]=", 89*x^2+32*x+38, " degree q=", 2
"f=", 94*x+94
"g=", 16*x^3+56*x^2+94*x
                                   5       4       3       2
        "Deterministic_NDSA:m: ", x  + 87 x  + 48 x  + 62 x  + 81 x + 29

"Deterministic_NDSA: u: ", 54*x^4+102*x^3+44*x^2+3*x+83
          "---------------------------------------------------------"

                                   "In MQRFR"

"q[", 1, "]=", 2*x+87, " degree q=", 1
"f=", 54*x^4+102*x^3+44*x^2+3*x+83
"g=", 1
"q[", 3, "]=", 19*x^2+40*x+64, " degree q=", 2
"f=", 94*x+71
"g=", 94*x^2+8*x+44
  "______________________________________________________________________________"

                      "MRFI sigma_[", 4, "] = ", [16, 81]

                            "In Deterministic_NDSA"

               "Deterministic_NDSA: alpha: ", [2, 3, 4, 5, 6, 7]

"Deterministic_NDSA:Psi_alpha: ",

    [[2, 90], [3, 97], [4, 104], [5, 4], [6, 11], [7, 18]]

"Deterministic_NDSA: Y: ",

    [[34, 20], [43, 100], [85, 49], [104, 46], [45, 80], [104, 49]]

                              6       5       4       3       2
   "Deterministic_NDSA:m: ", x  + 80 x  + 81 x  + 47 x  + 75 x  + 104 x + 11

"Deterministic_NDSA: u: ", 11*x^5+52*x^4+99*x^3+32*x^2+78*x+21
          "---------------------------------------------------------"

                                   "In MQRFR"

"q[", 1, "]=", 39*x+105, " degree q=", 1
"f=", 11*x^5+52*x^4+99*x^3+32*x^2+78*x+21
"g=", 1
"q[", 4, "]=", 13*x^2+91*x+45, " degree q=", 2
"f=", 80*x+80
"g=", 25*x^3+88*x^2+80*x
                                   5       4       3       2
        "Deterministic_NDSA:m: ", x  + 87 x  + 48 x  + 62 x  + 81 x + 29

"Deterministic_NDSA: u: ", x^4+105*x^3+58*x^2+84*x+48
          "---------------------------------------------------------"

                                   "In MQRFR"

"q[", 1, "]=", x+89, " degree q=", 1
"f=", x^4+105*x^3+58*x^2+84*x+48
"g=", 1
"q[", 3, "]=", 83*x^2+50*x+69, " degree q=", 2
"f=", 100*x+32
"g=", 100*x^2+31*x+106
  "______________________________________________________________________________"

                      "MRFI sigma_[", 5, "] = ", [32, 29]

                            "In Deterministic_NDSA"

               "Deterministic_NDSA: alpha: ", [2, 3, 4, 5, 6, 7]

"Deterministic_NDSA:Psi_alpha: ",

    [[2, 33], [3, 40], [4, 47], [5, 54], [6, 61], [7, 68]]

"Deterministic_NDSA: Y: ",

    [[12, 42], [102, 41], [16, 11], [89, 61], [26, 99], [49, 104]]

                              6       5       4       3       2
   "Deterministic_NDSA:m: ", x  + 80 x  + 81 x  + 47 x  + 75 x  + 104 x + 11

"Deterministic_NDSA: u: ", 80*x^5+59*x^4+19*x^3+18*x^2+3*x+23
          "---------------------------------------------------------"

                                   "In MQRFR"

"q[", 1, "]=", 103*x+20, " degree q=", 1
"f=", 80*x^5+59*x^4+19*x^3+18*x^2+3*x+23
"g=", 1
"q[", 4, "]=", 87*x^2+64*x+49, " degree q=", 2
"f=", 21*x+21
"g=", 40*x^3+78*x^2+21*x
                                   5       4       3       2
        "Deterministic_NDSA:m: ", x  + 87 x  + 48 x  + 62 x  + 81 x + 29

"Deterministic_NDSA: u: ", 5*x^4+73*x^3+70*x^2+77*x+14
          "---------------------------------------------------------"

                                   "In MQRFR"

"q[", 1, "]=", 43*x+53, " degree q=", 1
"f=", 5*x^4+73*x^3+70*x^2+77*x+14
"g=", 1
"q[", 3, "]=", 75*x^2+106*x+62, " degree q=", 2
"f=", 63*x+55
"g=", 63*x^2+64*x+9
  "______________________________________________________________________________"

                      "MRFI sigma_[", 6, "] = ", [64, 87]

                            "In Deterministic_NDSA"

               "Deterministic_NDSA: alpha: ", [2, 3, 4, 5, 6, 7]

"Deterministic_NDSA:Psi_alpha: ",

    [[2, 81], [3, 88], [4, 95], [5, 102], [6, 2], [7, 9]]

"Deterministic_NDSA: Y: ",

    [[22, 32], [28, 8], [29, 105], [16, 27], [92, 33], [86, 67]]

                              6       5       4       3       2
   "Deterministic_NDSA:m: ", x  + 80 x  + 81 x  + 47 x  + 75 x  + 104 x + 11

"Deterministic_NDSA: u: ", 102*x^5+69*x^4+32*x^3+68*x^2+15*x+18
          "---------------------------------------------------------"

                                   "In MQRFR"

"q[", 1, "]=", 64*x+54, " degree q=", 1
"f=", 102*x^5+69*x^4+32*x^3+68*x^2+15*x+18
"g=", 1
"q[", 4, "]=", 13*x^2+16*x+85, " degree q=", 2
"f=", 80*x+80
"g=", 25*x^3+10*x^2+80*x
                                   5       4       3       2
        "Deterministic_NDSA:m: ", x  + 87 x  + 48 x  + 62 x  + 81 x + 29

"Deterministic_NDSA: u: ", 90*x^4+46*x^3+100*x^2+28*x+15
          "---------------------------------------------------------"

                                   "In MQRFR"

"q[", 1, "]=", 44*x+51, " degree q=", 1
"f=", 90*x^4+46*x^3+100*x^2+28*x+15
"g=", 1
"q[", 3, "]=", 30*x^2+7*x+89, " degree q=", 2
"f=", 5*x+93
"g=", 5*x^2+2*x+16
  "______________________________________________________________________________"

                      "MRFI sigma_[", 7, "] = ", [21, 47]

                            "In Deterministic_NDSA"

               "Deterministic_NDSA: alpha: ", [2, 3, 4, 5, 6, 7]

"Deterministic_NDSA:Psi_alpha: ",

    [[2, 21], [3, 28], [4, 35], [5, 42], [6, 49], [7, 56]]

"Deterministic_NDSA: Y: ",

    [[61, 100], [81, 62], [26, 1], [21, 22], [61, 64], [17, 29]]

                              6       5       4       3       2
   "Deterministic_NDSA:m: ", x  + 80 x  + 81 x  + 47 x  + 75 x  + 104 x + 11

"Deterministic_NDSA: u: ", 91*x^5+56*x^4+42*x^3+48*x^2+2*x+1
          "---------------------------------------------------------"

                                   "In MQRFR"

"q[", 1, "]=", 20*x+65, " degree q=", 1
"f=", 91*x^5+56*x^4+42*x^3+48*x^2+2*x+1
"g=", 1
"q[", 4, "]=", 89*x^2+94*x+56, " degree q=", 2
"f=", 13*x+13
"g=", 91*x^3+91*x^2+13*x
                                   5       4       3       2
        "Deterministic_NDSA:m: ", x  + 87 x  + 48 x  + 62 x  + 81 x + 29

"Deterministic_NDSA: u: ", 2*x^4+43*x^3+80*x^2+6*x+34
          "---------------------------------------------------------"

                                   "In MQRFR"

"q[", 1, "]=", 54*x+6, " degree q=", 1
"f=", 2*x^4+43*x^3+80*x^2+6*x+34
"g=", 1
"q[", 3, "]=", 12*x^2+74*x+75, " degree q=", 2
"f=", 104*x+28
"g=", 104*x^2+104*x+76
  "______________________________________________________________________________"

"MRFI num_eval: ",

    [[92, 31, 16, 93, 33, 20, 101, 49], [107, 92, 47, 19, 42, 4, 104, 83]]

"MRFI den_eval: ",

    [[92, 2, 67, 34, 45, 28, 61, 75], [92, 1, 97, 31, 63, 41, 16, 80]]

  "______________________________________________________________________________"

                       "numerator_done: ", [false, false]

                      "denominator_done: ", [false, false]

  "--------------------------------------------------------------------------------"

"MRFI k=", 1
"MRFI numerator_done[", 1, "]=", false
  "______________________________________________________________________________"

                                   "In BMEA"

"v=", [92, 31, 16, 93, 33, 20, 101, 49]
"MRFI lambda_num: ", Z^2+104*Z+2
"MRFI terms_num: ", 2
"MRFI R_num: ", [[1, 1], [2, 1]]
  "--------------------------------------------------------------------------------"

                 "Checking termination condition for numerator"

                          "nops(R_num[", 1, "]): ", 2

                      "R_num[", 1, "]: ", [[1, 1], [2, 1]]

                           "terms_num[", 1, "]: ", 2

                                    "T: ", 4

                 "MRFI Numerator component ", 1, " recovered!"

  "____________________________________________________________________________"

"MRFI k=", 2
"MRFI numerator_done[", 2, "]=", false
  "______________________________________________________________________________"

                                   "In BMEA"

"v=", [107, 92, 47, 19, 42, 4, 104, 83]
"MRFI lambda_num: ", Z^2+103*Z+3
"MRFI terms_num: ", 2
"MRFI R_num: ", [[1, 1], [3, 1]]
  "--------------------------------------------------------------------------------"

                 "Checking termination condition for numerator"

                          "nops(R_num[", 2, "]): ", 2

                      "R_num[", 2, "]: ", [[1, 1], [3, 1]]

                           "terms_num[", 2, "]: ", 2

                                    "T: ", 4

                 "MRFI Numerator component ", 2, " recovered!"

  "____________________________________________________________________________"

  "===================================================================================="

                       "MFRI processing denominators now"

"MRFI den_eval: ", [[92, 2, 67, 34, 45, 28, 61, 75], [92, 1, 97, 31, 63, 41, 16
, 80]]
"MRFI common_den_flag: ", false
  "______________________________________________________________________________"

                                   "In BMEA"

"v=", [92, 2, 67, 34, 45, 28, 61, 75]
  "--------------------------------------------------------------------------------"

                "Checking termination condition for denominator"

                           "terms_den[", 1, "]: ", 2

                                    "T: ", 4

                "MRFI Denominator component ", 1, " recovered!"

  "______________________________________________________________________________"

                                   "In BMEA"

"v=", [92, 1, 97, 31, 63, 41, 16, 80]
  "--------------------------------------------------------------------------------"

                "Checking termination condition for denominator"

                           "terms_den[", 2, "]: ", 2

                                    "T: ", 4

                "MRFI Denominator component ", 2, " recovered!"

  "--------------------------------------------------------------------------------"

                       "BMEA done status: ", [true, true]

                           "All done status: ", true

"MRFI All components recovered!"
                   "MRFI Roots_num_eval: ", [[1, 2], [1, 3]]

                   "MRFI Roots_den_eval: ", [[2, 12], [1, 6]]

                     "-----------------------------------"

                            "In generate_monomials"

                               "roots_=", [1, 2]

                               "temp= ", [1, y1]

                       "MRFI num_mono[", 1, "]:", [1, y1]

                                                2
                  "MRFI lambda_num[", 1, "]:", Z  + 104 Z + 2

"In Zippel_Transpose_Vandermonde_solver"
                     "-----------------------------------"

                            "In generate_monomials"

                               "roots_=", [1, 3]

                               "temp= ", [1, y2]

                       "MRFI num_mono[", 2, "]:", [1, y2]

                                                2
                  "MRFI lambda_num[", 2, "]:", Z  + 103 Z + 3

"In Zippel_Transpose_Vandermonde_solver"
                            "In generate_monomials"

                               "roots_=", [2, 12]

"In Zippel_Transpose_Vandermonde_solver"
                            "In generate_monomials"

                               "roots_=", [1, 6]

"In Zippel_Transpose_Vandermonde_solver"
                        "In construct_final_polynomial"

                        "In construct_final_polynomial"

                        "In construct_final_polynomial"

                        "In construct_final_polynomial"

"MRFI ========================================"
"MRFI RECOVERY COMPLETE"
"MRFI ========================================"
> Ratrecon_num:=table():
> Ratrecon_den:=table():
> Final_rat_poly:=table():
> for i from 1 to num_eqn do 
>     Ratrecon_num[i]:=iratrecon(Num[i],p):
>     print("numerator = ",Ratrecon_num[i]):
>     Ratrecon_den[i]:=iratrecon(Den[i],p):
>     print("denominator = ",Ratrecon_den[i]):
> end do:
                             "numerator = ", 1 + y1

                                             2
                         "denominator = ", y1  y2 + y1

                            "numerator = ", -1 + y2

                          "denominator = ", y1 y2 + 1



> print("======================================================"):
            "======================================================"

> print("Displaying the results"):
                            "Displaying the results"


> if(num_eqn >1)then 
>     og_soln:=get_eqn(Sys,Vars):
>     # og_unordered_soln:=convert(og_soln,list):
>     # og_soln:=reording(og_unordered_soln,nops(Sys)):
>     # fin_rat_recon:=Vector(convert(Rat_recon,list)):
>     og_soln:=convert(og_soln,list):
>     for i from 1 to num_eqn do 
>         print("x",i,"="):
>         Final_rat_poly[i]:=Ratrecon_num[i]/Ratrecon_den[i]:
>         lprint("Rat_recon= ",Final_rat_poly[i]):
>         lprint("original_soln =",op(2,og_soln[i])):
>         #print("f",i,"/g",i,"-","ff",i,"/gg",i,"=",simplify(Rat_recon[i]-op(2,og_soln[i])));
>         printf("f%d/g%d-ff%d/gg%d = %a\n",i,i,i,i,simplify(Final_rat_poly[i]-op(2,og_soln[i])));
>     end do:
>     elif num_eqn =1 then 
>         Final_rat_poly[1]:=Ratrecon_num[1]/Ratrecon_den[1]:
>         lprint("Rat_recon= ",Final_rat_poly[1]):
>         lprint("Original polynomial =",F/G):
>         printf("f1/g1 - F/G = %a\n",simplify(Final_rat_poly[1]-F/G));
> end if:
                                  "in get_eqn"

                                  "x", 1, "="

"Rat_recon= ", (1+y1)/(y1^2*y2+y1)
"original_soln =", (1+y1)/y1/(y1*y2+1)
f1/g1-ff1/gg1 = 0
                                  "x", 2, "="

"Rat_recon= ", (-1+y2)/(y1*y2+1)
"original_soln =", (-1+y2)/(y1*y2+1)
f2/g2-ff2/gg2 = 0

> print("======================================================"):
            "======================================================"

> print("Total number of lines generated in get_point_on_affine_line:", num_lines):
       "Total number of lines generated in get_point_on_affine_line:", 9

> lprint("Total Black Box Calls:", counter):
"Total Black Box Calls:", 54

# param1:=6:
# param2:=2:
# debug_point:=B([param1,param2],p):
# print("B([",param1,",",param2,"],p) = ",debug_point):
# For param1=6, param2=2 we get "B([", 6, ",", 2, "],p) = ", [1366580498, 1561806289]
#  where igcd(1366580498, 1561806289)=1
# For rvr_bspline/main.mpl
# "B([", 6, ",", 2, "],p) = ", [976128928, 195225786]
# igcd(976128928, 195225786) = 2 and the reduced form is [488064464, 97612893]
# LEts try 111
# param1:=111:
# param2:=111:
# debug_point:=B([param1,param2],p):
# print("B([",param1,",",param2,"],p) = ",debug_point):
# print("igcd(debug_point[1], debug_point[2]) = ", igcd(debug_point[1], debug_point[2]) mod p):
> quit
memory used=9.1MB, alloc=41.3MB, time=0.07
