pixi global uninstall pypgcf
rm -rf ~/.pixi/envs/pypgcf
pixi clean --build
pixi clean cache --conda --yes
pixi clean cache --build --yes

pixi build -c -o ../
pixi global install --path $CONDA_FILE --channel conda-forge --channel bioconda --environment pypgc
