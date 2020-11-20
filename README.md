``` r
library(lemon) 
library(knitr)
```

I. Getting Started
==================

Installation
------------

To use `abmR`, you must first install it from Github using `devtools`
and load the library:

``` r
#devtools: install_github("bgoch5/abmr")
library(abmr,quietly=TRUE,warn.conflicts=FALSE)
```

While this package is still in development, it will be updated
frequently, so please be sure to re-install frequently.

For `abmR` to work, you will need updated versions of the packages
`geosphere`, `ggplot2`, `raster`, `rgdal`, `rnaturalearth`,
`rnaturalearthdata`, `sp`, `swfscMisc`, and `table1`.

Requirements for Environmental Data Rasters and Species objects
---------------------------------------------------------------

The modeling functions discussed below require several input objects to
work, which we will discuss here.

1.  A shapefile that represents the boundaries for your movements.

Here, we use North America, but whatever region you are interested in
modeling movement over.

``` r
library(rgdal) #For reading in shapefile
setwd("C:\\Users\\BGOCHANOUR\\Documents\\GitHub\\move-model")
NOAM=readOGR(layer="NOAM",dsn=".")
#> OGR data source with driver: ESRI Shapefile 
#> Source: "C:\Users\BGOCHANOUR\Documents\GitHub\move-model", layer: "NOAM"
#> with 36 features
#> It has 18 fields
plot(NOAM,xlim=c(-130,-60),ylim=c(18,75))
```

![](C:/Users/BGOCHANOUR/AppData/Local/Temp/RtmpEF3fhv/preview-33e4354d3935.dir/Vignette_files/figure-markdown_github/unnamed-chunk-3-1.png)

1.  A environmental raster stack

We are currently working on how to attach an example raster to the
package, but you can read yours in like this. Here, we use NDVI, but you
may use a raster stack for any environmental variable of interest. Here,
we plot the first four layers of our stack.

``` r
library(raster)
ndvi_raster=stack("C:\\Users\\BGOCHANOUR\\Documents\\GitHub\\move-model\\NDVI_2013.gri")
plot(ndvi_raster[[1:4]])
```

![](C:/Users/BGOCHANOUR/AppData/Local/Temp/RtmpEF3fhv/preview-33e4354d3935.dir/Vignette_files/figure-markdown_github/unnamed-chunk-4-1.png)

1.  A species object to describe the species that you will model.

This object contains information on origin longitude (x) and latitude
(y) as well as two morphological parameters, `p1` and `p2`. Use ?species
for more details.

``` r
# Define the species class
species <- setClass("species", slots=c(x="numeric", y="numeric",
p1="numeric", p1mean="numeric", p1sd="numeric",p1sign="character",
p2="numeric", p2mean="numeric", p2sd="numeric",p2sign="character"))

# Define your species object. Imagine p1 is wing chord and p2 is mass.
# We are simulating larger than average birds which increases movement speed.
# For more details, use ?species

my_species=species(x=-90,y=45,p1=15,p1mean=10,p1sd=2,p1sign="Pos",
                   p2=10,p2mean=8,p2sd=1,p2sign="Pos")
```

``` r
# Small fig.width
include_graphics("OCWA.jpg")
```

<img src="OCWA.jpg" width="50%" height="50%" />

II. Running Agent-based models
==============================

`abmR` has two functions for modeling animal migration: `moveSIM` and
`energySIM`. Both `moveSIM` and `energySIM` will output two things: (1)
a dataframe of results and (2) a summary table of input parameters. For
more details on how each function works, please see the documentation by
using ?moveSIM or ?energySIM.

`moveSIM`
---------

``` r
moveSIM_results=moveSIM(replicates=1,days=27,env_rast=ndvi_raster, search_radius=375,
  sigma=.4, dest_x=-100, dest_y=25, mot_x=1, mot_y=1,modeled_species=my_species,
  my_shapefile=NOAM,optimum=.5,direction="S",write_results=TRUE,single_rast=FALSE)
#> [1] "Agent died"
head(moveSIM_results$results,30)
```

|       lon|     lat|  day| agent\_id |    distance|
|---------:|-------:|----:|:----------|-----------:|
|   -90.000|  45.000|    1| Agent\_1  |          NA|
|   -89.575|  43.375|    2| Agent\_1  |   184.04681|
|   -90.225|  33.225|    3| Agent\_1  |  1131.30804|
|   -86.325|  31.525|    4| Agent\_1  |   412.57499|
|   -92.925|  29.825|    5| Agent\_1  |   659.50058|
|   -94.175|  30.625|    6| Agent\_1  |   149.62015|
|   -89.725|  29.675|    7| Agent\_1  |   441.17993|
|   -92.025|  30.175|    8| Agent\_1  |   228.76925|
|   -97.275|  29.275|    9| Agent\_1  |   517.26683|
|  -101.325|  23.175|   10| Agent\_1  |   790.18724|
|  -100.175|  26.275|   11| Agent\_1  |   364.14736|
|  -100.225|  24.425|   12| Agent\_1  |   206.00246|
|   -99.275|  24.825|   13| Agent\_1  |   105.94668|
|  -100.675|  26.125|   14| Agent\_1  |   201.83143|
|   -99.425|  24.625|   15| Agent\_1  |   209.01459|
|  -100.625|  25.925|   16| Agent\_1  |   188.50170|
|  -100.025|  24.875|   17| Agent\_1  |   131.53852|
|   -99.675|  24.625|   18| Agent\_1  |    45.01609|
|  -100.975|  24.975|   19| Agent\_1  |   137.02447|
|        NA|      NA|   20| Agent\_1  |          NA|
|        NA|      NA|   21| Agent\_1  |          NA|
|        NA|      NA|   22| Agent\_1  |          NA|
|        NA|      NA|   23| Agent\_1  |          NA|
|        NA|      NA|   24| Agent\_1  |          NA|
|        NA|      NA|   25| Agent\_1  |          NA|
|        NA|      NA|   26| Agent\_1  |          NA|
|        NA|      NA|   27| Agent\_1  |          NA|

``` r
moveSIM_results$run_params
```

<table style="width:100%;">
<colgroup>
<col style="width: 7%" />
<col style="width: 3%" />
<col style="width: 7%" />
<col style="width: 9%" />
<col style="width: 4%" />
<col style="width: 4%" />
<col style="width: 4%" />
<col style="width: 4%" />
<col style="width: 4%" />
<col style="width: 11%" />
<col style="width: 5%" />
<col style="width: 6%" />
<col style="width: 9%" />
<col style="width: 8%" />
<col style="width: 8%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: right;">replicates</th>
<th style="text-align: right;">days</th>
<th style="text-align: left;">env_raster</th>
<th style="text-align: right;">search_radius</th>
<th style="text-align: right;">sigma</th>
<th style="text-align: right;">dest_x</th>
<th style="text-align: right;">dest_y</th>
<th style="text-align: right;">mot_x</th>
<th style="text-align: right;">mot_y</th>
<th style="text-align: left;">modeled_species</th>
<th style="text-align: right;">optimum</th>
<th style="text-align: left;">direction</th>
<th style="text-align: left;">write_results</th>
<th style="text-align: left;">single_rast</th>
<th style="text-align: right;">missing_pct</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: right;">1</td>
<td style="text-align: right;">27</td>
<td style="text-align: left;">env_raster</td>
<td style="text-align: right;">375</td>
<td style="text-align: right;">0.4</td>
<td style="text-align: right;">-100</td>
<td style="text-align: right;">25</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: left;">my_species</td>
<td style="text-align: right;">0.5</td>
<td style="text-align: left;">S</td>
<td style="text-align: left;">TRUE</td>
<td style="text-align: left;">FALSE</td>
<td style="text-align: right;">29.62963</td>
</tr>
</tbody>
</table>

`energySIM`
-----------

``` r
energySIM_results=energySIM(replicates=1,days=27,env_rast=ndvi_raster, search_radius=375,
                  sigma=.4, dest_x=-100, dest_y=25, mot_x=1, mot_y=1,
                  modeled_species=my_species, my_shapefile=NOAM,
                  optimum_lo=.4,optimum_hi=.6,init_energy=100,
                  direction="S",write_results=TRUE,single_rast=FALSE)
#> [1] "Edge Case 2"
#> [1] "Edge Case 3"
head(energySIM_results$results,30)
```

|      lon|     lat|  energy|  day| agent\_id |   distance|  delta\_energy|
|--------:|-------:|-------:|----:|:----------|----------:|--------------:|
|  -90.000|  45.000|     100|    1| Agent\_1  |         NA|             NA|
|  -81.425|  43.775|      85|    2| Agent\_1  |  695.29728|            -15|
|  -78.725|  38.425|      65|    3| Agent\_1  |  637.07405|            -20|
|  -83.225|  37.875|      55|    4| Agent\_1  |  398.62257|            -10|
|  -81.525|  34.725|      65|    5| Agent\_1  |  382.36789|             10|
|  -86.825|  31.825|      45|    6| Agent\_1  |  589.36666|            -20|
|  -89.225|  30.625|      55|    7| Agent\_1  |  264.64030|             10|
|  -84.775|  30.375|      40|    8| Agent\_1  |  427.70478|            -15|
|  -85.225|  29.925|      55|    9| Agent\_1  |   66.22466|             15|
|  -80.825|  28.225|      55|   10| Agent\_1  |  467.99621|              0|
|  -81.125|  25.425|      45|   11| Agent\_1  |  313.11555|            -10|
|       NA|      NA|      NA|   12| Agent\_1  |         NA|             NA|
|       NA|      NA|      NA|   13| Agent\_1  |         NA|             NA|
|       NA|      NA|      NA|   14| Agent\_1  |         NA|             NA|
|       NA|      NA|      NA|   15| Agent\_1  |         NA|             NA|
|       NA|      NA|      NA|   16| Agent\_1  |         NA|             NA|
|       NA|      NA|      NA|   17| Agent\_1  |         NA|             NA|
|       NA|      NA|      NA|   18| Agent\_1  |         NA|             NA|
|       NA|      NA|      NA|   19| Agent\_1  |         NA|             NA|
|       NA|      NA|      NA|   20| Agent\_1  |         NA|             NA|
|       NA|      NA|      NA|   21| Agent\_1  |         NA|             NA|
|       NA|      NA|      NA|   22| Agent\_1  |         NA|             NA|
|       NA|      NA|      NA|   23| Agent\_1  |         NA|             NA|
|       NA|      NA|      NA|   24| Agent\_1  |         NA|             NA|
|       NA|      NA|      NA|   25| Agent\_1  |         NA|             NA|
|       NA|      NA|      NA|   26| Agent\_1  |         NA|             NA|
|       NA|      NA|      NA|   27| Agent\_1  |         NA|             NA|

``` r
energySIM_results$run_params
```

<table style="width:100%;">
<colgroup>
<col style="width: 6%" />
<col style="width: 2%" />
<col style="width: 6%" />
<col style="width: 8%" />
<col style="width: 3%" />
<col style="width: 4%" />
<col style="width: 4%" />
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 9%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 5%" />
<col style="width: 8%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 8%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: right;">replicates</th>
<th style="text-align: right;">days</th>
<th style="text-align: left;">env_raster</th>
<th style="text-align: right;">search_radius</th>
<th style="text-align: right;">sigma</th>
<th style="text-align: right;">dest_x</th>
<th style="text-align: right;">dest_y</th>
<th style="text-align: right;">mot_x</th>
<th style="text-align: right;">mot_y</th>
<th style="text-align: left;">modeled_species</th>
<th style="text-align: right;">optimum_lo</th>
<th style="text-align: right;">optimum_hi</th>
<th style="text-align: left;">direction</th>
<th style="text-align: left;">write_results</th>
<th style="text-align: left;">single_rast</th>
<th style="text-align: right;">missing_pct</th>
<th style="text-align: right;">mortality_pct</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: right;">1</td>
<td style="text-align: right;">27</td>
<td style="text-align: left;">env_raster</td>
<td style="text-align: right;">375</td>
<td style="text-align: right;">0.4</td>
<td style="text-align: right;">-100</td>
<td style="text-align: right;">25</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: left;">my_species</td>
<td style="text-align: right;">0.4</td>
<td style="text-align: right;">0.6</td>
<td style="text-align: left;">S</td>
<td style="text-align: left;">TRUE</td>
<td style="text-align: left;">FALSE</td>
<td style="text-align: right;">59.25926</td>
<td style="text-align: right;">0</td>
</tr>
</tbody>
</table>

III. Visualizing Results
========================

The output from `moveSIM` and `energySIM` are simple dataframes, so one
may easily visualize their results using custom codes. However, we have
provided the following functions that include a set of default plots
that you may find useful.

`moveVIZ`
---------

One may want a basic plot of each agent’s movement, as compared to a
straight-line trajectory between the start and target endpoint. This can
be obtained by

``` r
moveVIZ(moveSIM_results,title="MoveSIM results")
```

![](C:/Users/BGOCHANOUR/AppData/Local/Temp/RtmpEF3fhv/preview-33e4354d3935.dir/Vignette_files/figure-markdown_github/unnamed-chunk-9-1.png)

`energyVIZ`
-----------

The `energyVIZ` function can produce three different types of output, so
you must specify the `type` parameter. For a similar basic plot to the
above, use `type` = “plot”

``` r
energyVIZ(energySIM_results,title="EnergySIM results",type="plot")
```

![](C:/Users/BGOCHANOUR/AppData/Local/Temp/RtmpEF3fhv/preview-33e4354d3935.dir/Vignette_files/figure-markdown_github/unnamed-chunk-10-1.png)

ADD GRADIENT/CONTOUR PLOT SOON

One may also want to create a table to summarize energy, day, and
distance, stratified by energy gain/loss amount. This can be obtained by
setting the `type` argument to “table”.

``` r
energyVIZ(energySIM_results,title="EnergySIM results",type="table")
#> [1] "<table class=\"Rtable1\">\n<thead>\n<tr>\n<th class='rowlabel firstrow lastrow'></th>\n<th class='firstrow lastrow'><span class='stratlabel'>-20<br><span class='stratn'>(N=2)</span></span></th>\n<th class='firstrow lastrow'><span class='stratlabel'>-15<br><span class='stratn'>(N=2)</span></span></th>\n<th class='firstrow lastrow'><span class='stratlabel'>-10<br><span class='stratn'>(N=2)</span></span></th>\n<th class='firstrow lastrow'><span class='stratlabel'>0<br><span class='stratn'>(N=1)</span></span></th>\n<th class='firstrow lastrow'><span class='stratlabel'>10<br><span class='stratn'>(N=2)</span></span></th>\n<th class='firstrow lastrow'><span class='stratlabel'>15<br><span class='stratn'>(N=1)</span></span></th>\n<th class='firstrow lastrow'><span class='stratlabel'>Overall<br><span class='stratn'>(N=27)</span></span></th>\n</tr>\n</thead>\n<tbody>\n<tr>\n<td class='rowlabel firstrow'><span class='varlabel'>energy</span></td>\n<td class='firstrow'></td>\n<td class='firstrow'></td>\n<td class='firstrow'></td>\n<td class='firstrow'></td>\n<td class='firstrow'></td>\n<td class='firstrow'></td>\n<td class='firstrow'></td>\n</tr>\n<tr>\n<td class='rowlabel'>Mean (SD)</td>\n<td>55.0 (14.1)</td>\n<td>62.5 (31.8)</td>\n<td>50.0 (7.07)</td>\n<td>55.0 (NA)</td>\n<td>60.0 (7.07)</td>\n<td>55.0 (NA)</td>\n<td>60.5 (18.0)</td>\n</tr>\n<tr>\n<td class='rowlabel'>Median [Min, Max]</td>\n<td>55.0 [45.0, 65.0]</td>\n<td>62.5 [40.0, 85.0]</td>\n<td>50.0 [45.0, 55.0]</td>\n<td>55.0 [55.0, 55.0]</td>\n<td>60.0 [55.0, 65.0]</td>\n<td>55.0 [55.0, 55.0]</td>\n<td>55.0 [40.0, 100]</td>\n</tr>\n<tr>\n<td class='rowlabel lastrow'>Missing</td>\n<td class='lastrow'>0 (0%)</td>\n<td class='lastrow'>0 (0%)</td>\n<td class='lastrow'>0 (0%)</td>\n<td class='lastrow'>0 (0%)</td>\n<td class='lastrow'>0 (0%)</td>\n<td class='lastrow'>0 (0%)</td>\n<td class='lastrow'>16 (59.3%)</td>\n</tr>\n<tr>\n<td class='rowlabel firstrow'><span class='varlabel'>day</span></td>\n<td class='firstrow'></td>\n<td class='firstrow'></td>\n<td class='firstrow'></td>\n<td class='firstrow'></td>\n<td class='firstrow'></td>\n<td class='firstrow'></td>\n<td class='firstrow'></td>\n</tr>\n<tr>\n<td class='rowlabel'>Mean (SD)</td>\n<td>4.50 (2.12)</td>\n<td>5.00 (4.24)</td>\n<td>7.50 (4.95)</td>\n<td>10.0 (NA)</td>\n<td>6.00 (1.41)</td>\n<td>9.00 (NA)</td>\n<td>14.0 (7.94)</td>\n</tr>\n<tr>\n<td class='rowlabel lastrow'>Median [Min, Max]</td>\n<td class='lastrow'>4.50 [3.00, 6.00]</td>\n<td class='lastrow'>5.00 [2.00, 8.00]</td>\n<td class='lastrow'>7.50 [4.00, 11.0]</td>\n<td class='lastrow'>10.0 [10.0, 10.0]</td>\n<td class='lastrow'>6.00 [5.00, 7.00]</td>\n<td class='lastrow'>9.00 [9.00, 9.00]</td>\n<td class='lastrow'>14.0 [1.00, 27.0]</td>\n</tr>\n<tr>\n<td class='rowlabel firstrow'><span class='varlabel'>distance</span></td>\n<td class='firstrow'></td>\n<td class='firstrow'></td>\n<td class='firstrow'></td>\n<td class='firstrow'></td>\n<td class='firstrow'></td>\n<td class='firstrow'></td>\n<td class='firstrow'></td>\n</tr>\n<tr>\n<td class='rowlabel'>Mean (SD)</td>\n<td>613 (33.7)</td>\n<td>562 (189)</td>\n<td>356 (60.5)</td>\n<td>468 (NA)</td>\n<td>324 (83.2)</td>\n<td>66.2 (NA)</td>\n<td>424 (187)</td>\n</tr>\n<tr>\n<td class='rowlabel'>Median [Min, Max]</td>\n<td>613 [589, 637]</td>\n<td>562 [428, 695]</td>\n<td>356 [313, 399]</td>\n<td>468 [468, 468]</td>\n<td>324 [265, 382]</td>\n<td>66.2 [66.2, 66.2]</td>\n<td>413 [66.2, 695]</td>\n</tr>\n<tr>\n<td class='rowlabel lastrow'>Missing</td>\n<td class='lastrow'>0 (0%)</td>\n<td class='lastrow'>0 (0%)</td>\n<td class='lastrow'>0 (0%)</td>\n<td class='lastrow'>0 (0%)</td>\n<td class='lastrow'>0 (0%)</td>\n<td class='lastrow'>0 (0%)</td>\n<td class='lastrow'>17 (63.0%)</td>\n</tr>\n</tbody>\n</table>\n"
```

IV. Getting Help
================

Need help? Have suggestions to make `abmR` better? If so, please open an
issue or pull request on Github
(<a href="https://github.com/bgoch5/abmr" class="uri">https://github.com/bgoch5/abmr</a>),
or drop me an email at ben.gochanour@ou.edu.



abmR: An R Package for Agent-based Model Analysis of Large-scale Movements Across Taxa
======================================================================================

#### Benjamin Gochanour, Javi Fernandez Lopez, Andrea Contina

#### 2020-10-26

Installation
------------

To use `abmR`, you must first install it from Github using `devtools`
and load the library:

``` {.r}
#devtools: install_github("bgoch5/abmr")
library(abmr,quietly=TRUE,warn.conflicts=FALSE)
```

For the full `abmR` vignette, see the following website:
https://www.bengochanour.com/vignetteupdates
