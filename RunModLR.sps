* Encoding: UTF-8.
* Encoding: UTF-8.

* Make sure the Active Dataset shows your data file (.sav).
* Steps:
* First. Open ModLR.sps, select all, and Run Selection (click the triangle green button in the menu above the syntax).
*         You have just invoked the macro in the SPSS memory, which you can call during your SPSS session.
* Second. Highlight the code below, and Run Selection (click the triangle green button above).


MLR iv = q8_6_speed
/covs=
/mod = gencat
/dv = q5_recom
/tVar =  0
/tLine =  0
/tiv =  0
/tcat =  1
/hci= 0
/ramsey= 0
/quadratic= 0
/multim= 0.

