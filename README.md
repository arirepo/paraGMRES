# Important Update not Covered in the Report!

When running the code, you need to provide an extra argument for the number of restarts! So, for example, a correct command to run the dense 6x6 matrix with 2 processes and one restart would be:

$> mpirun -np 2 ./gmres ../matrices/dense6x6.mtx void 6 0 1

Thanks to Sean Smith sesm630_AT_gmail.com for providing the feedback.

# paraGMRES
Massively Scalable Parallel GMRES C-code for Sparse System of Equations
![alt tag](https://raw.github.com/arrgasm/paraGMRES/master/banner.png)
