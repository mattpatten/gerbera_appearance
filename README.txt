This is the whole pipeline of cleaning and standardising flower images, identifying flower properties, loading in behavioural 
ratings and performing a regression analysis to identify the flower properties that explain aesthetic preferences.

The main script to cover all of this is run_pipeline.m and just requires a source directory containing the flower images.
Key scripts are in the base directory, more specific and helper scripts are located in appropriately-named folders for the
parts of the pipeline the scripts are aiding.


Briefly, the pipeline performs the following:
(i)   Rename the image to a standardized format.
(ii)  Finding the flower face in the image, and generate a mask for regions outside of this area.
(iii) Identify flower centre, pad/crop surround to make image a perfect square and rescale image to specific dimensions (e.g., 300x300).
(iv)  Perform transformations of images, as required (L*a*b*/grey colourspaces, applying mask, finding pixels along mask border).
(v)   Perform Kovesi's phase congruency algorithm to identify stimulus edges and high frequency areas.
(vi)  Perform k-means clustering on the L*a*b* colour channels to identify eccentricities with rapid colour changes.
(vii) Segment flower into disk florets, trans florets and ray florets based on the above two algorithms. 
(viii)Allow manual corrections to flower segmentations.
(ix)  Generate values for ~45 flower attributes (e.g., width of ray florets), displaying representative images and histograms.
(x)   Uses the values across all of these properties to group the flower set into a specific number of clusters.
(xi)  Performs a regression analysis on these properties in order to explain behavioural ratings measured elsewhere (see flwrpoll)
      to identify the attributes that contribute to human preference.

Created by Matt Patten on 27 Nov 2018
