# DissNet

DissNet is an open-source Python package for exploring, dissecting, analyzing
and visualizing human structural connectome data. It includes modules for data
visualization, Graph Theory measures extraction, standard connectivity measures
estimation. All the measures can be also computed for every single network.

<p align="center">
    <b>New Visualization Tool</b>
</p>  
<p align="center">
    <img src="https://github.com/Davi1990/DissNet/blob/main/docs/video_unscreen.gif" width="380"/>
</p>

<p align="center">
    <b>Network-level connectivity measures</b>
</p>
<p align="center">
    <img src="https://github.com/Davi1990/DissNet/blob/main/docs/network.png" width="380"/> <img src="https://github.com/Davi1990/DissNet/blob/main/docs/spider_plot.png" width="380"/>
</p>

<p align="center">
    <b>Network-level vulnerability against target or random attack</b>
</p>
<p align="center">
    <img src="https://github.com/Davi1990/DissNet/blob/main/docs/resilience_net.png" width="380"/> <img src="https://github.com/Davi1990/DissNet/blob/main/docs/vulnerability_net.png" width="380"/>
</p>


# Install
To install the latest stable version of DissNet, you can use setup.py in a terminal:

```bash
    cd DissNet/
    python setup.py install --user
```

To test that all has been installed correctly type in a terminal:

```bash
    cd DissNet/
    python test.py
```


# Dependencies
- Python>=3.5
- Pysurfer
- Pandas
- MNE
- Mayavi
- networkx
