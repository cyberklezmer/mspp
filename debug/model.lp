\ENCODING=ISO-8859-1
\Problem name: IloCplex

Minimize
 obj: x0@0#0 + u@0#1 + 0.5 theta@0m0#2 + 0.5 theta@0m1#3
Subject To
 c1:  x0@0#0 + x0@0m0#4 >= 0
 c2:  x0@0#0 + x0@0m0#4 <= 1
 c3:  10.5 u@0#1 + theta@0m0#2 + theta@0m1#3 - 8.4 x0@0m0#4 - 10.5 u@0m0#5
      - 5.25 theta@0m0m0#6 - 5.25 theta@0m0m1#7 - 8.4 x0@0m1#8 - 10.5 u@0m1#9
      - 5.25 theta@0m1m0#10 - 5.25 theta@0m1m1#11 >= 0
 c4:  0.5 u@0#1 + theta@0m0#2 + theta@0m1#3 - 0.4 x0@0m0#4 - 0.5 u@0m0#5
      - 0.25 theta@0m0m0#6 - 0.25 theta@0m0m1#7 - 0.4 x0@0m1#8 - 0.5 u@0m1#9
      - 0.25 theta@0m1m0#10 - 0.25 theta@0m1m1#11 >= 0
 c5:  x0@0#0 + x0@0m0#4 + x0@0m0m0#12 >= 0
 c6:  x0@0#0 + x0@0m0#4 + x0@0m0m0#12 <= 1
 c7:  10.5 u@0m0#5 + theta@0m0m0#6 + theta@0m0m1#7 + 21 u@0m1#9
      + theta@0m1m0#10 + theta@0m1m1#11 - 6.72 x0@0m0m0#12 - 10.5 u@0m0m0#13
      - 5.25 theta@0m0m0m0#14 - 5.25 theta@0m0m0m1#15 - 6.72 x0@0m0m1#16
      - 10.5 u@0m0m1#17 - 5.25 theta@0m0m1m0#18 - 5.25 theta@0m0m1m1#19
      - 6.72 x0@0m1m0#20 - 10.5 u@0m1m0#21 - 5.25 theta@0m1m0m0#22
      - 5.25 theta@0m1m0m1#23 - 6.72 x0@0m1m1#24 - 10.5 u@0m1m1#25
      - 5.25 theta@0m1m1m0#26 - 5.25 theta@0m1m1m1#27 >= 0
 c8:  0.5 u@0m0#5 + theta@0m0m0#6 + theta@0m0m1#7 + u@0m1#9 + theta@0m1m0#10
      + theta@0m1m1#11 - 0.32 x0@0m0m0#12 - 0.5 u@0m0m0#13
      - 0.25 theta@0m0m0m0#14 - 0.25 theta@0m0m0m1#15 - 0.32 x0@0m0m1#16
      - 0.5 u@0m0m1#17 - 0.25 theta@0m0m1m0#18 - 0.25 theta@0m0m1m1#19
      - 0.32 x0@0m1m0#20 - 0.5 u@0m1m0#21 - 0.25 theta@0m1m0m0#22
      - 0.25 theta@0m1m0m1#23 - 0.32 x0@0m1m1#24 - 0.5 u@0m1m1#25
      - 0.25 theta@0m1m1m0#26 - 0.25 theta@0m1m1m1#27 >= 0
 c9:  x0@0#0 + x0@0m0#4 + x0@0m0m0#12 + x0@0m0m0m0#28  = 1
 c10: 10.5 u@0m0m0#13 + theta@0m0m0m0#14 - 5.376 x0@0m0m0m0#28 >= 0
 c11: 0.5 u@0m0m0#13 + theta@0m0m0m0#14 - 0.256 x0@0m0m0m0#28 >= 0
 c12: x0@0#0 + x0@0m0#4 + x0@0m0m0#12 + x0@0m0m0m1#29  = 1
 c13: 10.5 u@0m0m0#13 + theta@0m0m0m1#15 - 7.392 x0@0m0m0m1#29 >= 0
 c14: 0.5 u@0m0m0#13 + theta@0m0m0m1#15 - 0.352 x0@0m0m0m1#29 >= 0
 c15: x0@0#0 + x0@0m0#4 + x0@0m0m1#16 >= 0
 c16: x0@0#0 + x0@0m0#4 + x0@0m0m1#16 <= 1
 c17: 10.5 u@0m0#5 + theta@0m0m1#7 + 21 u@0m1#9 + theta@0m1m0#10
      + theta@0m1m1#11 - 9.24 x0@0m0m1#16 - 10.5 u@0m0m1#17
      - 5.25 theta@0m0m1m0#18 - 5.25 theta@0m0m1m1#19 - 9.24 x0@0m1m0#20
      - 10.5 u@0m1m0#21 - 5.25 theta@0m1m0m0#22 - 5.25 theta@0m1m0m1#23
      - 9.24 x0@0m1m1#24 - 10.5 u@0m1m1#25 - 5.25 theta@0m1m1m0#26
      - 5.25 theta@0m1m1m1#27 >= 0
 c18: 0.5 u@0m0#5 + theta@0m0m1#7 + u@0m1#9 + theta@0m1m0#10 + theta@0m1m1#11
      - 0.44 x0@0m0m1#16 - 0.5 u@0m0m1#17 - 0.25 theta@0m0m1m0#18
      - 0.25 theta@0m0m1m1#19 - 0.44 x0@0m1m0#20 - 0.5 u@0m1m0#21
      - 0.25 theta@0m1m0m0#22 - 0.25 theta@0m1m0m1#23 - 0.44 x0@0m1m1#24
      - 0.5 u@0m1m1#25 - 0.25 theta@0m1m1m0#26 - 0.25 theta@0m1m1m1#27 >= 0
 c19: x0@0#0 + x0@0m0#4 + x0@0m0m1#16 + x0@0m0m1m0#30  = 1
 c20: 10.5 u@0m0m1#17 + theta@0m0m1m0#18 - 7.392 x0@0m0m1m0#30 >= 0
 c21: 0.5 u@0m0m1#17 + theta@0m0m1m0#18 - 0.352 x0@0m0m1m0#30 >= 0
 c22: x0@0#0 + x0@0m0#4 + x0@0m0m1#16 + x0@0m0m1m1#31  = 1
 c23: 10.5 u@0m0m1#17 + theta@0m0m1m1#19 - 10.164 x0@0m0m1m1#31 >= 0
 c24: 0.5 u@0m0m1#17 + theta@0m0m1m1#19 - 0.484 x0@0m0m1m1#31 >= 0
 c25: x0@0#0 + x0@0m1#8 >= 0
 c26: x0@0#0 + x0@0m1#8 <= 1
 c27: 10.5 u@0#1 + theta@0m1#3 - 11.55 x0@0m1#8 - 10.5 u@0m1#9
      - 5.25 theta@0m1m0#10 - 5.25 theta@0m1m1#11 >= 0
 c28: 0.5 u@0#1 + theta@0m1#3 - 0.55 x0@0m1#8 - 0.5 u@0m1#9
      - 0.25 theta@0m1m0#10 - 0.25 theta@0m1m1#11 >= 0
 c29: x0@0#0 + x0@0m1#8 + x0@0m1m0#20 >= 0
 c30: x0@0#0 + x0@0m1#8 + x0@0m1m0#20 <= 1
 c31: 10.5 u@0m1#9 + theta@0m1m0#10 + theta@0m1m1#11 - 9.24 x0@0m1m0#20
      - 10.5 u@0m1m0#21 - 5.25 theta@0m1m0m0#22 - 5.25 theta@0m1m0m1#23
      - 9.24 x0@0m1m1#24 - 10.5 u@0m1m1#25 - 5.25 theta@0m1m1m0#26
      - 5.25 theta@0m1m1m1#27 >= 0
 c32: 0.5 u@0m1#9 + theta@0m1m0#10 + theta@0m1m1#11 - 0.44 x0@0m1m0#20
      - 0.5 u@0m1m0#21 - 0.25 theta@0m1m0m0#22 - 0.25 theta@0m1m0m1#23
      - 0.44 x0@0m1m1#24 - 0.5 u@0m1m1#25 - 0.25 theta@0m1m1m0#26
      - 0.25 theta@0m1m1m1#27 >= 0
 c33: x0@0#0 + x0@0m1#8 + x0@0m1m0#20 + x0@0m1m0m0#32  = 1
 c34: 10.5 u@0m1m0#21 + theta@0m1m0m0#22 - 7.392 x0@0m1m0m0#32 >= 0
 c35: 0.5 u@0m1m0#21 + theta@0m1m0m0#22 - 0.352 x0@0m1m0m0#32 >= 0
 c36: x0@0#0 + x0@0m1#8 + x0@0m1m0#20 + x0@0m1m0m1#33  = 1
 c37: 10.5 u@0m1m0#21 + theta@0m1m0m1#23 - 10.164 x0@0m1m0m1#33 >= 0
 c38: 0.5 u@0m1m0#21 + theta@0m1m0m1#23 - 0.484 x0@0m1m0m1#33 >= 0
 c39: x0@0#0 + x0@0m1#8 + x0@0m1m1#24 >= 0
 c40: x0@0#0 + x0@0m1#8 + x0@0m1m1#24 <= 1
 c41: 10.5 u@0m1#9 + theta@0m1m1#11 - 12.705 x0@0m1m1#24 - 10.5 u@0m1m1#25
      - 5.25 theta@0m1m1m0#26 - 5.25 theta@0m1m1m1#27 >= 0
 c42: 0.5 u@0m1#9 + theta@0m1m1#11 - 0.605 x0@0m1m1#24 - 0.5 u@0m1m1#25
      - 0.25 theta@0m1m1m0#26 - 0.25 theta@0m1m1m1#27 >= 0
 c43: x0@0#0 + x0@0m1#8 + x0@0m1m1#24 + x0@0m1m1m0#34  = 1
 c44: 10.5 u@0m1m1#25 + theta@0m1m1m0#26 - 10.164 x0@0m1m1m0#34 >= 0
 c45: 0.5 u@0m1m1#25 + theta@0m1m1m0#26 - 0.484 x0@0m1m1m0#34 >= 0
 c46: x0@0#0 + x0@0m1#8 + x0@0m1m1#24 + x0@0m1m1m1#35  = 1
 c47: 10.5 u@0m1m1#25 + theta@0m1m1m1#27 - 13.9755 x0@0m1m1m1#35 >= 0
 c48: 0.5 u@0m1m1#25 + theta@0m1m1m1#27 - 0.6655 x0@0m1m1m1#35 >= 0
Bounds
 0 <= x0@0#0 <= 1
      u@0#1 Free
      theta@0m0#2 >= -1000000
      theta@0m1#3 >= -1000000
      u@0m0#5 Free
      theta@0m0m0#6 Free
      theta@0m0m1#7 Free
      u@0m1#9 Free
      theta@0m1m0#10 Free
      theta@0m1m1#11 Free
      u@0m0m0#13 Free
      theta@0m0m0m0#14 Free
      theta@0m0m0m1#15 Free
      u@0m0m1#17 Free
      theta@0m0m1m0#18 Free
      theta@0m0m1m1#19 Free
      u@0m1m0#21 Free
      theta@0m1m0m0#22 Free
      theta@0m1m0m1#23 Free
      u@0m1m1#25 Free
      theta@0m1m1m0#26 Free
      theta@0m1m1m1#27 Free
End
