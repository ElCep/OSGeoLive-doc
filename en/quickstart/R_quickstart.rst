:Author: Barry Rowlingson, Etienne Delay
:Version: osgeo-live11.0
:License: Creative Commons Attribution-ShareAlike 3.0 Unported  (CC BY-SA 3.0)

.. image:: /images/project_logos/logo-R.png
  :alt: project logo
  :align: right

********************************************************************************
R Quickstart
********************************************************************************

R is a free software environment for statistical computing and graphics.
A greate value of R is his strong commmunity providing usfull documentation.

We warmly encourage you to take a look on the `R cookbook <http://www.cookbook-r.com/>`_ for basic use.

In this quickstart we will focus on somes of R spatial possibility :

  * load some data from a shapefile or raster
  * do a coordinate transformation
  * basic uses of ggplot2 with spatial data
  *

Start R
================================================================================

Either:

  * Choose ``R Statistics`` from the Spatial Tools section of the  menu - a terminal window appears running R.
  * Enter ``R`` at a command-line shell prompt. R will start in that terminal.

Don't fear the command line - it is a source of great power. Using the up and down arrows to recall commands so you can edit mistakes will help greatly.
Hit CTRL-C if get stuck and you should get the prompt back.

R through the `The Comprehensive R Archive Network <https://cran.r-project.org/>`_ provide a huge number of packages helping users to work. In this quikstart we will propose to use several packages not basically installed.

If a missing `library()` is called in a script, you can download and install it with `install.packages("name_library")`.

Quit R
-------------------------------------------------------------------------------

Almost everything in R is a function, including the function for quitting. Type
``q()`` and hit return. If you just type ``q`` you'll see the source code for the ``q`` function.

R will ask you if you want to save your workspace as an R data image file. When you
start R again from a directory with a ``.RData`` file it will restore all its
data from there.


Beginning R
-------------------------------------------------------------------------------

R is essentially a command-line program, although graphical user
interfaces are available. You type a line of code at the prompt,
press return, and the R interpreter evaluates it and prints the
result.

You can start with simple arithmetic

::

   > 3*2
   [1] 6

   > 1 + 2 * 3 / 4
   [1] 2.5

   > sqrt(2)
   [1] 1.414214

   > pi * exp(-1)
   [1] 1.155727


And so on. A full range of arithmetic, trigonometric, and statistical
functions are built in, and thousands more are available from
packages in the `CRAN <http://cran.r-project.org/>`_ archive.

The main prompt in R is ``>``, but there is also the continuation prompt, ``+``, which appears if R expects more input to make a valid expression.
You'll see this if you forget a closing bracket or parenthesis.

::

   > sqrt(
   + 2
   + )
   [1] 1.414214




Loading vector data
-------------------------------------------------------------------------------

There are many packages for spatial data manipulation and statistics. Some
are included here, and some can be downloaded from CRAN.

Here we will load two shapefiles - the country boundaries and populated places
from the Natural Earth data. We use two add-on packages to get the spatial
functionality:

::

	> library(sp)
	> library(maptools)

	> countries = readShapeSpatial("/usr/local/share/data/natural_earth2/ne_10m_admin_0_countries.shp")
	> places = readShapeSpatial("/usr/local/share/data/natural_earth2/ne_10m_populated_places.shp")
	> plot(countries)

This gives us a simple map of the world:

.. image:: /images/screenshots/1024x768/r_plot1.png

If _maptools_ packages provide useful function, you can also use _rgdal_ package and `writeOGR() <https://www.rdocumentation.org/packages/rgdal/versions/1.2-8/topics/writeOGR>`_

When an OGR dataset is read into R in this way we get back an object that
behaves in many ways like a data frame. We can use the ``admin``
field to subset the world data and just get the UK:

::

	> uk = countries[countries$admin == "United Kingdom",]
	> plot(uk); axis(1); axis(2)

.. image:: /images/screenshots/1024x768/r_plot2.png

This looks a bit squashed to anyone who lives here, since we are more familiar with
a coordinate system centred at our latitude. Currently the object doesn't have a
coordinate system assigned to it - we can check this with some more functions:

::

	> proj4string(uk)
	[1] NA

``NA`` is a missing data marker. We need to assign a CRS to the object before we can
transform it with the spTransform function from the rgdal package. We transform
to EPSG:27700 which is the Ordnance Survey of Great Britain grid system:

::

	> proj4string(uk) = CRS("+init=epsg:4326")
	> library(rgdal)
	> ukos = spTransform(uk, CRS("+init=epsg:27700"))
	> proj4string(ukos)
	[1] " +init=epsg:27700 +proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs
	+towgs84=446.448,-125.157,542.060,0.1502,0.2470,0.8421,-20.4894"

	> plot(ukos); axis(1); axis(2)

This plots the base map of the transformed data. Now we want to add some points from the
populated place data set. Again we subset the points we want and transform them to
Ordnance Survey Grid Reference coordinates:

::

	> ukpop = places[places$ADM0NAME == "United Kingdom",]
	> proj4string(ukpop) = CRS("+init=epsg:4326")
	> ukpop = spTransform(ukpop, CRS("+init=epsg:27700"))

We add these points to the base map, scaling their size by scaled square root of the
population (because that makes a symbol with area proportional to population), set the
colour to red and the plotting character to a solid blob:
::

	> points(ukpop, cex=sqrt(ukpop$POP_MAX/1000000), col="red", pch=19)
	> title("UK Population centre sizes")

and our final image appears:

.. image:: /images/screenshots/1024x768/r_plot3.png

Load and map raster data
-------------------------------------------------------------------------------
First we download raster from `MODIS <https://en.wikipedia.org/wiki/Moderate-resolution_imaging_spectroradiometer>`_ data using R.

::

    library(MODIS)
    #modis option modified on .MODIS_Opts.R file on ~
    #donwload and tile preparation of all NDVI data:

    dates = as.POSIXct( as.Date(c("01/04/2016","30/06/2016"),format = "%d/%m/%Y") )
    dates2 = transDate(dates[1],dates[2]) # Transform input dates from before
    # The MODIS package allows you to select tiles interactively.  We, however, select them manually here
    runGdal("MOD13Q1",begin=dates2$beginDOY,end=dates2$endDOY,
            tileH=16:16,tileV=7,outProj="+init=epsg:32628")

This will download of each images the hdf and proceed to the post-processing. At the end you will find 12 raster for each date :

* 16_days_NDVI
* 16_days_EVI
* 16_days_VI_Quality
* 16_days_red_reflectance
* 16_days_NIR_reflectance
* 16_days_blue_reflectance
* 16_days_MIR_reflectance
* 16_days_view_zenith_angle
* 16_days_sun_zenith_angle
* 16_days_relative_azimuth_angle
* 16_days_composite_day_of_the_year
* 16_days_pixel_reliability.tif



::

    library(raster) ## raster manipulation and stat
    library(sp) ## spatial function (proj, resolution etc.)
    library(doParallel) ## parallelisation some action
    library(rgdal) ## another way to load shape file
    library(magrittr) ## use pipe in R



    #create stack from tif NDVI files
    file.l = list.files(path = "~/MODIS_ARC/PROCESSED/MOD13Q1.006_20170826140906/", pattern = "250m_16_days_NDVI")


    my.stack = stack(paste0(path.v,l.files)) ##create raster stack of files in a directory
    # load borders from previous file
    countries = readShapeSpatial("/usr/local/share/data/natural_earth2/ne_10m_admin_0_countries.shp")
    border = countries[countries@data$ADMIN == "Senegal",]
    border = spTransform(border, "+init=epsg:32628")

Now you can take a look at your data with `plot()`. We have used the stack to create a lite. To plot just one of these NDVI layer we use the double square parenthesis. The `add = TRUE` allow us to add some data to the first plot.

::

    plot(my.stack[[1]])
    plot(border, add = TRUE)

.. _NDVI_border:
.. image:: /images/screenshots/1024x768/r-ndvi_sng.png

Now we want to crop all of our raster within the same time with the `border` layer. We can easily use `doParallel` package to parallels this task.

::

    ## We use doParallel and magrittr packages to pipe different actions (crop and mask)
    registerDoParallel(6) #we will use 6 parallel thread
    result = foreach(i = 1:dim(my.stack)[3],.packages='raster',.inorder=T) %dopar% {
      my.stack[[i]] %>%
        crop(border) %>%
        mask(border)
    }
    endCluster()
    ndvi.stack = stack(result)

You have see on the previous map (NDVI_border_.) than unity are not what we are waiting for NDVI values (something between -1 and 1). So let's resale our data !

::

    ndvi.stack = ndvi.stack*0.0001 #rescaling of MODIS data
    ndvi.stack[ndvi.stack ==-0.3]=NA #Fill value(-0,3) in NA
    ndvi.stack[ndvi.stack < (-0.2)]=NA # as valide range is -0.2 -1 , all values smaller than -0,2 are masked out
    names(ndvi.stack) = seq.POSIXt(from = ISOdate(2016,4,1), by = "16 day", length.out = 6)

And finaly let's plot it !

::

    my_palette = colorRampPalette(c("red", "yellow", "lightgreen")) #Create a color palette for our values
    ## Plot X maps in the same layout
    spplot(ndvi.stack, layout=c(2, 3),
        col.regions = my_palette(16))

.. image:: /images/screenshots/1024x768/r_sng_senegal_crop.png

Plots with ggplot2
================================================================================


Dynamic map with leaflet
================================================================================


Vignettes
================================================================================

In the past the documentation for R packages tended to be tersely-written help pages for each function. Now package authors are encouraged to write a 'vignette' as a friendly introduction to the package.

If you just run the ``vignette()`` function with no arguments you will get the list of those vignettes on your system. Try ``vignette("intro_sp")`` for a slightly technical introduction to the R spatial package.

The ``vignette("gstat")`` gives a tutorial in the use of that package for spatial interpolation including Kriging.

.. comment: doesn't work// or ``vignette("shapefiles")`` for explanations of using shapefiles in R.

Further Reading
================================================================================

* For general information about R, try the official `Introduction to R <http://cran.r-project.org/doc/manuals/R-intro.html>`_ or any of the documentation from the main `R Project <http://www.r-project.org/>`_ page.
* You can also take time with the `R cookbook <http://www.cookbook-r.com/>`_. A great place to learn basis and go deeper into mechanics of some indispensable packages like *ggplot2* or *ddplyr*.
* For more information on spatial aspects of R, the best place to start is probably the `R Spatial Task View <http://cran.r-project.org/web/views/Spatial.html>`_
* RStoolbox for remote sensing in R have a `nice blog <http://bleutner.github.io/RStoolbox/>`_

R and space in books
===============================================================================
* Wegmann, Martin, Benjamin Leutner, et Stefan Dech, éd. 2016. _Remote Sensing and GIS for Ecologists: Using Open Source Software_. Pelagic Publishing.

* Blangiardo, Marta, et Michela Cameletti. 2015. _Spatial and Spatio-temporal Bayesian Models with R - INLA_. 1ʳᵉ éd. Wiley.

* Brunsdon, Chris, et Lex Comber. 2015. _An Introduction to R for Spatial Analysis & Mapping_. Los Angeles: SAGE Publications Ltd.
