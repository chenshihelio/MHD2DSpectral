install:
    git clone https://github.com/chenshihelio/MHD2DSpectral.git
    micromamba env create --file environment.yml
    micromamba activate laps
    make

clean:
    rm *.o *.mod