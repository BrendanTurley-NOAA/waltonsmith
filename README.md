# waltonsmith
Processing R/V Walton Smith data

Update 2022-06-01: Bottom data interpolation now covaries with bathymetry

The interpolation is done using the [fields](https://cran.r-project.org/web/packages/fields/index.html) R package which has nice kriging methods and handy visualization functions. There is an excellent [vignette](https://github.com/NCAR/fields/blob/master/fieldsVignette.pdf) that explains the main functions and features.

Underdevelopment:
- interpolation that takes into consideration coastlines
- add city names to maps to make them recognizable
- add polygons to blank out inland waters

---

# Latest cruise plot

![alt text](https://github.com/imaginaryfish/waltonsmith/blob/main/figures/latest_underway.png "latest underway data")

![alt text](https://github.com/imaginaryfish/waltonsmith/blob/main/figures/latest_bottom.png "latest bottom data")