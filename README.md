# KSMM
Kernel Support Matrix Machine

```ksmm``` is an R package. ```ksmm``` provides functions to fit support matrix machine, and it also provides datas that used to fit.

## 1. INSTALLATION

```ksmm``` is not submitted to CRAN. Therefore, the ```ksmm``` can not be installed through install.packages("ksmm") in R prompt.
Instead, ```ksmm``` can be installed through our GitHub.
Please see below to install in R.

(1) From GitHub
```{r}
> library(devtools)
> install_github("staycuriousn/ksmm")
```

## 2. USAGE NOTES

(1) Access to and use of real data in the manuscript

- Nottingham data : The nottingham data is the collection of images useful for conducting experiments in faces. This data can be available at [pics.stir.ac.uk](pics.stir.ac.uk)
```{r}
> library(ksmm)
> data(nottingham)
```
- Brainwave data : The EEG signal data from 10 college students while they watched MOOC video clips. Click the [link](https://www.kaggle.com/wanghaohan/confused-eef) for more details.
```{r}
> library(ksmm)
> data(braninwave)
```
- EEG Alcohol data : The multiple electrode time recordings of control and alcoholic subjects. Click the [link](http://kdd.ics.uci.edu/databases/eeg/eeg.data.html) for more details.
```{r}
> library(ksmm)
> data(EEG_alcohol)
```

(2) Description of R functions in ```ksmm``` 

- Descriptions of arguments in the functions in ```ksmm``` can be obtained by ```help()``` or ```?``` in R prompt, and documentation of ```ksmm```.   


(3) List of R functions in ```ksmm```

- ```ksmm``` : ```ksmm``` function is used to fit the type I and II ZILNBGM model, and ZILPGM model.

- ```kfold_ksmm``` : ```kfold_ksmm``` function tunes hyperparameters of statistical methods using a grid search over supplied parameter ranges for ksmm.


