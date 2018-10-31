# Q_truncate

The package was tested on TACC's Stampede2 using 

* gcc/7.1.0
* impi/18.0.2
* Boost/1.64 (MultiArray)
* Eigen3/3.3.5 (Included in the package)

To compile the code:

```
cd src
make QTRMPI=1
```

To run the program (Stampede2 SKX node):

```
module load gcc/7.1.0
module load impi/18.0.2
module load boost/1.64

export OMP_NUM_THREADS=1

time ibrun ./qtr_mpi ini.test
```

Example of input file (2-D Scattering):

```
# ini.test

[MAIN]
job=scatter2d
inFilename=ini.test
quiet=true
[SCATTERXD]
dimensions=2
isFullGrid=false
period=100
k=0.001
h1=0.05
h2=0.05
xi1=-4
xf1=8
xi2=-4
xf2=8
Tf=50
hb=1
m=1
TolH=1e-7
TolL=1e-6
TolHd=1e-3
TolLd=5e-3
ExReduce=0
lambda=0.11803
```




