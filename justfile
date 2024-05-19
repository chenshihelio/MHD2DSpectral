# set the default number of processors (this number should be 2^n)
np := "4"

default: 
    just --list

env:
    micromamba env create --file environment.yml

install:
    cd src && make install && make clean

[no-cd]
run np=np:
    mpirun -np {{np}} lasp
    mkdir -p diags
    mv *.dat diags/

[no-cd]
clean:
    find . -name 'log' -type f -delete
    find . -name '*.o' -type f -delete
    find . -name '*.mod' -type f -delete
    find . -name '*.dat' -type f -delete
    find . -name '.DS_Store' -type f -delete