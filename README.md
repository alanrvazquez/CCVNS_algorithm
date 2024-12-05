README file for MATLAB code to accompany Vazquez, A. R., Goos, P., and Schoen E. D. (2019). "Constructing 
Two-Level Designs by Concatenation of Strength-3 Orthogonal Arrays." Technometrics, Vol. 61, p.p. 219-232.

===CONTENTS=======================================

A. SCRIPT FILES
B. PRIMARY FUNCTIONS
C. UTILITY FUNCTIONS
D. DATA FILES
==================================================

A) SCRIPT FILES.
Use these to construct concatenated designs using the CC/VNS algorithm.

Examples_B4.m           Shows how to construct an concatenated designs that optimizes
                        the B4 value.
                 
Examples_F4.m           Shows how to construct an concatenated designs that sequentially
                        optimizes the F4 vector.


B) PRIMARY FUNCTIONS.
These are the functions one would interact with when optimizing and evaluating 
the concatenated designs. 

CCVNSCFV.m              Standard version of the CC/VNS algorithm for optimizing the 
                        F4 vector of concatenated designs of strength three.
                        
Parallel_CCVNSCFV.m     Parallel version of the CC/VNS algorithm for optimizing the 
                        F4 vector of concatenated designs of strength three. It 
                        distributes the iterations of the algorithm into multiple 
                        cores available in your computer. Useful for fast computations. 
                        We recommend this version.
                        
CCVNSBfour.m            Standard version of the CC/VNS algorithm for optimizing the 
                        B4 value of concatenated designs of strength three.
                        
Parallel_CCVNSBfour.m   Parallel version of the CC/VNS algorithm for optimizing the 
                        B4 value of concatenated designs of strength three. It 
                        distributes the iterations of the algorithm into multiple 
                        cores available in your computer. Useful for fast computations. 
                        
F4.m                    Computes the F4 vector and B4 value of a two-level orthogonal 
                        array.


C) UTILITY FUNCTIONS.
These are called by the primary functions for specific purposes.

CCAlgBfour.m            Column Change (CC) algorithm for optimizing the B4 value.

Bfour.m                 Computes the B4 value of two-level orthogonal arrays using 
                        Butler's moment matrix. For details see: Butler, N. A. (2003). 
                        Minimum aberration construction results for nonregular two-level 
                        factorial designs. Biometrika, 90:891-898.

bigMFfour.m             Computes the f(D) values of a design. For more information, see
                        Supplementary Section A.2.
                        
CCAlgCFV.m              Column Change (CC) algorithm for the optimizing the f(D) function 
                        based on the F4 vector. 
                        
ChangeColImpact.m       Calculates the F4 vector of a concatenated design constructed from a 
                        lower design that swaps columns 'ii' and 'jj'. Useful for improving
                        computing times when constructing and evaluating plans for 
                        the lower design. For more information, see Supplementary Section A.2

SignSwitchImpact.m      Calculates the J4-characteristics over 4-factor sets involving 
                        a specific column. Useful for improving the computing times when 
                        testing a sign switch in a column in the lower design. For more 
                        information, see Supplementary Section A.2.
                        
TwoFIMat.m              Computes the two-factor interaction matrix of a design.
                        This function is faster than the x2fx function in Matlab.

D) DATA FILES.

OAN32K16.mat            Complete catalog of strength-3 orthogonal arrays with 32 runs 
                        and 16 factors. Obtained from Schoen, E. D., Eendebak P. T., 
                        and Nguyen M. V. M. (2010). Complete enumeration of pure-level 
                        and mixed-level orthogonal arrays. Journal of Combinatorial 
                        Designs, 18:-123-140. For more details, please visit
                        http://www.pietereendebak.nl/oapackage/. 

