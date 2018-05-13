param num_matrices >= 1, integer; # Number of matrices in the data file to be read
param num_rows >= 1, integer;     # Number of rows
param num_cols >= 1, integer;     # Number of columns 

set MATS    := 1 .. num_matrices; # set of matrices
set ROWS    := 1 .. num_rows;	  # set of rows
set COLUMNS := 1 .. num_cols;	  # set of columns

param matrix_value {MATS, ROWS, COLUMNS} >= 0; # values for entries of each matrix

param tumor_value{ROWS,COLUMNS} >= 0;

param critical_value{ROWS,COLUMNS} >= 0;

set tumor_area := {k in ROWS, h in COLUMNS: tumor_value[k,h] > 0};
set critical_area := {k in ROWS, h in COLUMNS: critical_value[k,h] > 0};

var S {(k,h) in tumor_area} >= 0;
var T {(k,h) in critical_area} >= 0;

param tumor_min = 10;
param tumor_max = 70;
param critical_max = 2;

# Full Matrix of 1's
set N1 = 1 .. 60;
set N2 = 1 .. 80;
param full_matrix {i in N1,j in N2} >=0;
let {i in N1,j in N2} full_matrix[i,j] := 1;

# Non Critical Area
set NC_area := {k in ROWS, h in COLUMNS: full_matrix[k,h] > 0
						&& (k,h) not in tumor_area
						&& (k,h) not in critical_area
					  };

# Dosage scalar of each beam
var X {i in MATS} >= 0;
param q = 2;
param w = 1;

# Minimize total Dosage: Non Critical And Critical Dosage
minimize total_dosage : q* sum {(k,h) in critical_area} sum {i in MATS} X[i] * matrix_value[i,k,h] + w*sum {(k,h) in NC_area} sum {i in MATS} X[i] * matrix_value[i,k,h] + sum {(k,h) in critical_area} T[k,h]; 

# Constraints: Tumor dose upper and lower limit, Critical Dose relaxed
s.t. tumor_limit_min {(k,h) in tumor_area} : sum {i in MATS} X[i] * matrix_value[i,k,h] >= tumor_min;

s.t. tumor_limit_max {(k,h) in tumor_area} : sum {i in MATS} X[i] * matrix_value[i,k,h] <= tumor_max;

s.t. critical_limit {(k,h) in critical_area} : sum {i in MATS} X[i] * matrix_value[i,k,h] == critical_max + T[k,h];


