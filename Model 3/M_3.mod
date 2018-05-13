param num_matrices >= 1, integer; # Number of matrices in the data file to be read
param num_rows >= 1, integer;     # Number of rows
param num_cols >= 1, integer;     # Number of columns 

set MATS    := 1 .. num_matrices; # set of matrices
set ROWS    := 1 .. num_rows;	  # set of rows
set COLUMNS := 1 .. num_cols;	  # set of columns

# values for entries of each matrix
param target_value {MATS, ROWS, COLUMNS} >= 0; 
param tumor{ROWS,COLUMNS} >= 0;
param critical{ROWS,COLUMNS} >= 0;

set tumor_area := {m in ROWS, n in COLUMNS: tumor[m,n] > 0};
set critical_area := {m in ROWS, n in COLUMNS: critical[m,n] > 0};

var A {(m,n) in tumor_area} >= 0;
var B {(m,n) in critical_area} >= 0;

#given lower limit for tumor dosage
param tumor_lower = 10;

#given upper limit for critical dosage
param critical_upper = 2;

#penalising factor
param q = 5;

# The coordinates of min and max values
param minval_rows = min {(m,n) in critical_area } m-1; 
param minval_col = min {(m,n) in critical_area } n-1; 

param maxval_rows = max {(m,n) in critical_area } m+1; 
param maxval_col = max {(m,n) in critical_area } n+1; 

# Full Matrix 
set N1 = 1 .. 60;
set N2 = 1 .. 80;
param full_matrix {i in N1,j in N2} >=0;
let {i in N1,j in N2} full_matrix[i,j] := 1;

# Border area
set border_area := {m in ROWS, n in COLUMNS: full_matrix[m,n] > 0 
						&& n >= minval_col
						&& n <= maxval_col
						&& m >= minval_rows
						&& m <= maxval_rows
						&&(m,n) not in tumor_area
						&&(m,n) not in critical_area
					  };

# intensity multiplier of each beam
var X {i in MATS} >= 0;

# Minimize critical Dosage till the limit
minimize total : sum {(m,n) in critical_area} sum {i in MATS} X[i] * target_value[i,m,n] +  sum {(m,n) in critical_area} B[m,n] -
sum {(m,n) in tumor_area} sum {i in MATS} X[i] * target_value[i,m,n] + sum {(m,n) in tumor_area} A[m,n] + q*sum {(m,n) in border_area} sum {i in MATS} X[i] * target_value[i,m,n]; 

#Subject ot the following constraints
s.t. border {(m,n) in border_area}: sum {i in MATS} X[i] * target_value[i,m,n] <= critical_upper;

s.t. tumor_limit {(m,n) in tumor_area} : sum {i in MATS} X[i] * target_value[i,m,n] == tumor_lower - A[m,n];

s.t. critical_limit {(m,n) in critical_area} : sum {i in MATS} X[i] * target_value[i,m,n] == critical_upper + B[m,n];




