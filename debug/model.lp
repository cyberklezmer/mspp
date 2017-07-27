\ENCODING=ISO-8859-1
\Problem name: IloCplex

Minimize
 obj: x0@0#0 + u@0#1 + 0.5 theta@0m0#2 + 0.5 u@0m0#3 + 0.25 theta@0m0m0#4
      + 0.25 u@0m0m0#5 + 0.125 theta@0m0m0m0#6 + 0.125 theta@0m0m0m1#7
      + 0.25 theta@0m0m1#8 + 0.25 u@0m0m1#9 + 0.125 theta@0m0m1m0#10
      + 0.125 theta@0m0m1m1#11 + 0.5 theta@0m1#12 + 0.5 u@0m1#13
      + 0.25 theta@0m1m0#14 + 0.25 u@0m1m0#15 + 0.125 theta@0m1m0m0#16
      + 0.125 theta@0m1m0m1#17 + 0.25 theta@0m1m1#18 + 0.25 u@0m1m1#19
      + 0.125 theta@0m1m1m0#20 + 0.125 theta@0m1m1m1#21
Subject To
 c1:  x0@0#0 + x0@0m0#22 >= 0
 c2:  x0@0#0 + x0@0m0#22 <= 1
 c3:  10.5 u@0#1 + theta@0m0#2 - 8.4 x0@0m0#22 >= 0
 c4:  0.5 u@0#1 + theta@0m0#2 - 0.4 x0@0m0#22 >= 0
 c5:  x0@0#0 + x0@0m0#22 + x0@0m0m0#23 >= 0
 c6:  x0@0#0 + x0@0m0#22 + x0@0m0m0#23 <= 1
 c7:  10.5 u@0m0#3 + theta@0m0m0#4 - 6.72 x0@0m0m0#23 >= 0
 c8:  0.5 u@0m0#3 + theta@0m0m0#4 - 0.32 x0@0m0m0#23 >= 0
 c9:  x0@0#0 + x0@0m0#22 + x0@0m0m0#23 + x0@0m0m0m0#24  = 1
 c10: 10.5 u@0m0m0#5 + theta@0m0m0m0#6 - 5.376 x0@0m0m0m0#24 >= 0
 c11: 0.5 u@0m0m0#5 + theta@0m0m0m0#6 - 0.256 x0@0m0m0m0#24 >= 0
 c12: x0@0#0 + x0@0m0#22 + x0@0m0m0#23 + x0@0m0m0m1#25  = 1
 c13: 10.5 u@0m0m0#5 + theta@0m0m0m1#7 - 7.392 x0@0m0m0m1#25 >= 0
 c14: 0.5 u@0m0m0#5 + theta@0m0m0m1#7 - 0.352 x0@0m0m0m1#25 >= 0
 c15: x0@0#0 + x0@0m0#22 + x0@0m0m1#26 >= 0
 c16: x0@0#0 + x0@0m0#22 + x0@0m0m1#26 <= 1
 c17: 10.5 u@0m0#3 + theta@0m0m1#8 - 9.24 x0@0m0m1#26 >= 0
 c18: 0.5 u@0m0#3 + theta@0m0m1#8 - 0.44 x0@0m0m1#26 >= 0
 c19: x0@0#0 + x0@0m0#22 + x0@0m0m1#26 + x0@0m0m1m0#27  = 1
 c20: 10.5 u@0m0m1#9 + theta@0m0m1m0#10 - 7.392 x0@0m0m1m0#27 >= 0
 c21: 0.5 u@0m0m1#9 + theta@0m0m1m0#10 - 0.352 x0@0m0m1m0#27 >= 0
 c22: x0@0#0 + x0@0m0#22 + x0@0m0m1#26 + x0@0m0m1m1#28  = 1
 c23: 10.5 u@0m0m1#9 + theta@0m0m1m1#11 - 10.164 x0@0m0m1m1#28 >= 0
 c24: 0.5 u@0m0m1#9 + theta@0m0m1m1#11 - 0.484 x0@0m0m1m1#28 >= 0
 c25: x0@0#0 + x0@0m1#29 >= 0
 c26: x0@0#0 + x0@0m1#29 <= 1
 c27: 10.5 u@0#1 + theta@0m1#12 - 11.55 x0@0m1#29 >= 0
 c28: 0.5 u@0#1 + theta@0m1#12 - 0.55 x0@0m1#29 >= 0
 c29: x0@0#0 + x0@0m1#29 + x0@0m1m0#30 >= 0
 c30: x0@0#0 + x0@0m1#29 + x0@0m1m0#30 <= 1
 c31: 10.5 u@0m1#13 + theta@0m1m0#14 - 9.24 x0@0m1m0#30 >= 0
 c32: 0.5 u@0m1#13 + theta@0m1m0#14 - 0.44 x0@0m1m0#30 >= 0
 c33: x0@0#0 + x0@0m1#29 + x0@0m1m0#30 + x0@0m1m0m0#31  = 1
 c34: 10.5 u@0m1m0#15 + theta@0m1m0m0#16 - 7.392 x0@0m1m0m0#31 >= 0
 c35: 0.5 u@0m1m0#15 + theta@0m1m0m0#16 - 0.352 x0@0m1m0m0#31 >= 0
 c36: x0@0#0 + x0@0m1#29 + x0@0m1m0#30 + x0@0m1m0m1#32  = 1
 c37: 10.5 u@0m1m0#15 + theta@0m1m0m1#17 - 10.164 x0@0m1m0m1#32 >= 0
 c38: 0.5 u@0m1m0#15 + theta@0m1m0m1#17 - 0.484 x0@0m1m0m1#32 >= 0
 c39: x0@0#0 + x0@0m1#29 + x0@0m1m1#33 >= 0
 c40: x0@0#0 + x0@0m1#29 + x0@0m1m1#33 <= 1
 c41: 10.5 u@0m1#13 + theta@0m1m1#18 - 12.705 x0@0m1m1#33 >= 0
 c42: 0.5 u@0m1#13 + theta@0m1m1#18 - 0.605 x0@0m1m1#33 >= 0
 c43: x0@0#0 + x0@0m1#29 + x0@0m1m1#33 + x0@0m1m1m0#34  = 1
 c44: 10.5 u@0m1m1#19 + theta@0m1m1m0#20 - 10.164 x0@0m1m1m0#34 >= 0
 c45: 0.5 u@0m1m1#19 + theta@0m1m1m0#20 - 0.484 x0@0m1m1m0#34 >= 0
 c46: x0@0#0 + x0@0m1#29 + x0@0m1m1#33 + x0@0m1m1m1#35  = 1
 c47: 10.5 u@0m1m1#19 + theta@0m1m1m1#21 - 13.9755 x0@0m1m1m1#35 >= 0
 c48: 0.5 u@0m1m1#19 + theta@0m1m1m1#21 - 0.6655 x0@0m1m1m1#35 >= 0
Bounds
 0 <= x0@0#0 <= 1
      u@0#1 Free
      theta@0m0#2 Free
      u@0m0#3 Free
      theta@0m0m0#4 Free
      u@0m0m0#5 Free
      theta@0m0m0m0#6 Free
      theta@0m0m0m1#7 Free
      theta@0m0m1#8 Free
      u@0m0m1#9 Free
      theta@0m0m1m0#10 Free
      theta@0m0m1m1#11 Free
      theta@0m1#12 Free
      u@0m1#13 Free
      theta@0m1m0#14 Free
      u@0m1m0#15 Free
      theta@0m1m0m0#16 Free
      theta@0m1m0m1#17 Free
      theta@0m1m1#18 Free
      u@0m1m1#19 Free
      theta@0m1m1m0#20 Free
      theta@0m1m1m1#21 Free
End
