# Use of machine learning methods forsuperpixel-based spatial classification and segmentation
Engineer's thesis code supplement
------

In this repistory one can find the implementation of ideas presented in my engineering thesis together with the data needed to recreate the results.

All scripts can be found in the scripts folder:
- *install_packages.R* - run the code in this script to install all the required packages
- *source.R* - script with functions used in other scripts, mostly technical
- *lrs_hrs.R* - script with the implementation of my method
- *smst.R* - script with the implementation of the [SMST method](https://ieeexplore.ieee.org/abstract/document/8038855)

In both *lrs_hrs.R* and *smst.R* the image to be processed needs to be chosen.
All values needed for both approaches are stored in the *choose_image* variable.

Choosing the image is done by supplying the image date or image number for row extraction from the choose_image variable:
```R

current_image = choose_image["image-date",]

#or 

current_image = choose_image[image_number,]

```

As an example:
```R

current_image = choose_image["2017-05-19",]

#or 

current_image = choose_image[1,]

```

