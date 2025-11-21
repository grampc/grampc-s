This directory contains the source files for the documentation. The documentation is build via sphinx. 

First, you have to install python. Then It's recommended to create a virtual environment via
```
python -m venv .venv
```
or one can use a context action from vscode.
Afterwards, install the requirements inside the requirements.txt file located in this folder with
```
pip install -r requirements.txt
```

To build the doc, type `make html` for html output or `make latexpdf` for pdf output. 