Need to install

https://github.com/MTG/sms-tools

This includes most of the course material as well.

To work with freesound, use https://github.com/MTG/freesound-python

------

To install and run ipython notebook on cloud9.  Note that cloud9 sessions often die, so you'll have to restart the notebook (so save frequently).

1. Download+install miniconda:
  * from http://jupyter.readthedocs.org/en/latest/install.html#new-to-python-and-jupyter
  * then https://www.continuum.io/downloads
  * then http://conda.pydata.org/miniconda.html (for miniconda)
  * then pick the appropriate link, for linux it's a bash script:
    ```bash
    $ curl https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh > x.sh
    $ chmod +x x.sh
    $ ./x.sh
    ```
    (there will be some prompts that you have to answer)
    At the end of the output, it will say where it installed miniconda, probably someplace like /home/ubuntu/miniconda2/bin/conda

1. Install jpyter:
 ```bash
/home/ubuntu/miniconda2/bin/conda install jupyter 
```
 
1. Install numpy + scipy (optional)
 ```bash
/home/ubuntu/miniconda2/bin/conda install numpy scipy
```

1. Run the notebook:

 ```bash
/home/ubuntu/miniconda2/bin/jupyter notebook --ip=0.0.0.0 --port=8080 --no-browser
```

1. Access the notebook at:

 ```bash
https://<workspace>-<username>.c9users.io/
```

 e.g.

 ```bash
https://asp-conanyuan.c9users.io/
```
