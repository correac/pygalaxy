pygalaxy
=======

[pygalaxy](https://github.com/correac/pygalaxy) is an analytic formalism
based on Correa et al. (2018a,b) to calculate the best fitting expression
of hot halo gas heating & cooling rates, gas mass histories, hot gas fractions.


Note that pygalaxy assumes halo virial mass (M200) is 200 times critical overdensity, and
concentration is the ratio of halo virial mass (R200) over scale radius (obtained from best-fit NFW profile)

Requirements
------------

Written in python, it uses routines in numpy and scipy to create a structured dataset. The package requires:

+ `python3.6` or above
+ see requirements.txt


Installing
----------

To get started using the package you need to set up a python virtual environment. The steps are as follows:

Clone pygalaxy
```
git clone https://github.com/correac/pygalaxy.git

cd pygalaxy

python3 -m venv pygalaxy_env
```

Now activate the virtual environment.

```
source pygalaxy_env/bin/activate
```

Update pip just in case
```
pip install pip --upgrade

pip install -r requirements.txt
```

How to use it
-------------

To run the script type
```bash
 python3 pygalaxy.py -z input_redshift \
                     --calculate-mcrit yes \
                     -o path_to_output_directory 
```

### Support and Contact

If you have trouble with pygalaxy or you have feature requests/suggestions please
open an issue at https://github.com/correac/pygalaxy/issues
