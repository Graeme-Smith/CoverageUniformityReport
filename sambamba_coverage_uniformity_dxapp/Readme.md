# sambamba_coverage_uniformity (DNAnexus Platform App)

## What does this app do?
This app uses the calculated Sambamba coverage bed files for x samples to visualise the variability in coverage at each genomic region across all these samples.  Results are available as an interactive plot, simplified pdf plot, and csv table.

## What are typical use cases for this app?

This app can be used during the development of new diagnositic tests to identify regions which have low coverage.

## What inputs are required for this app to run?
This app requires the following inputs:
 - Name of project containing the samabamba output folder coverage/raw_output/. 

## How does this app work?
At the core of the app is a dockerised R script which reads in all data from the folder coverage/raw_output/ in the project, cleans up the data, and outputs html, pdf, and csv reports.

## What does this app output?

This app outputs:
 - An inteactive HTML report with boxplots of every region sorted by coverage from low > high
 - A static plot showing the average for each region
 - A csv file with the region data in tabular form

## What are the limitations of this app?

 - The app currently does not work with files created with MokaCan - the prescence of any of these output files in the input folder will cause the script to stall as it tries to plot > 100,000 data points.
 - Currently outputs are not defined in the dxapp.json causing an error message when the app completes - the output folder is instead uploaded via a hardcoded line in the shell script and placed in coverage/uniformity_metrics/ within the named project.
 - Currently you cannot define input files - the script looks in a predefined folder for the input and this filepath is hard-coded into the script.

