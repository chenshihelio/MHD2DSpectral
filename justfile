default: 
    just --list

env:
    micromamba env create --file environment.yml

install:
    cd src && make install && make clean

run:
    mpirun -np 4 lasp
    mkdir -p {{invocation_directory()}}/diags
    mv *.dat {{invocation_directory()}}/diags/

clean:
    find . -name 'log' -type f -delete
    find . -name '*.o' -type f -delete
    find . -name '*.mod' -type f -delete
    find . -name '*.dat' -type f -delete
    find . -name '.DS_Store' -type f -delete