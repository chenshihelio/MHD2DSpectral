install:
    micromamba env create --file environment.yml
    micromamba activate laps
    make

clean:
    rm *.o *.mod