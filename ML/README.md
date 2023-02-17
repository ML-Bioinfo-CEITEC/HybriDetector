# CNN for smallRNA binding prediction

Convolutional Neural Network for prediction of smallRNA:target site binding based on given smallRNA and target sequence.


### Web application

The user-friendly web application for miRNA and tRNA target predictions https://ml-bioinfo-ceitec.github.io/HybriDetector/


## Installation of standalone version

Using Git:
```
git clone https://github.com/ML-Bioinfo-CEITEC/HybriDetector.git
```

### Prerequisites

It is implemented in python using Keras and Tensorflow backend.

Required:

* python, recommended version 3.7
    * Keras 2.7.0
    * tensorflow 2.7.0
    * pandas
    * numpy
    
    
### Installing

```
# create a virtual environment:

python -m venv env

# activate it and install the necessary libraries.

source env/bin/activate
cd ML
pip install -r requirements.txt
```

### Prediction

To make the predictions, choose one of the miRNA (Models/model_miRNA.h5), tRNA (Models/model_tRNA.h5) or yRNA (Models/model_yRNA.h5) models. Required input is a tsv file with multiple potential small RNA - target pairs consisting of first column containing small RNA sequence (idealy around 20 bp long) and second column containing exactly 50 bp long target sequence.

```
# if you are not actively sourcing from the previously created virtualenv:
source venv/bin/activate
# run the prediction
./predict.py --input <input_file> --output <output_file> --model <Models/model_miRNA.h5>
```

