# Seven Bridges Challenge

This project assumes that your machine had already *pip*, *virtualenv*. If you do not have them, you can check the links below for instructions:
* [pip](https://packaging.python.org/installing/#install-pip-setuptools-and-wheel)
* [virtualenv](https://python-guide-pt-br.readthedocs.io/en/latest/dev/virtualenvs/)

In order to start project you will do the following (Assuming you are on the right directory):

```
cd sbchallenge
virtualenv env
source env/bin/activate
pip install -r requirements.txt
chmod +x kmer_counter.py
./kmer_counter.py "filename" "k_mer size" "top count"
```
