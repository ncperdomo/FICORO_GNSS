# Filtering and Combining GNSS Velocity Fields in Python

## 1) Introduction

The Mediterranean and the Middle East regions provide a unique opportunity to investigate plate dynamics and large-scale crustal deformation and earthquake hazard given the complex tectonic interactions involving multiple plates and tectonic regimes. The expansion of regional GNSS networks and the availability of published velocities have facilitated a broader understanding of active tectonics in these regions. However, despite these advancements, few attempts have been made to integrate the available GNSS velocities, and no consensus exists in the geodetic community on the most effective way to filter and harmonise geodetic data sets. The primary objective of this project is to filter, combine and rotate GNSS velocity fields from column-formatted text files. These data are obtained from different geodetic studies of lithospheric deformation, with each file comprising 13 distinct columns, including 'Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj', 'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig', 'Stat'. The goal is to combine all the published GNSS velocity fields in the Euro-Mediterranean region, Asia Minor, and the Middle East, rotate them into different reference frames, filter outliers, and manage repeated stations that might differ in names or have coordinates varying by $\leq 0.01^\circ$ (1.11 km) across different studies. By providing insights into plate motions and crustal deformation, our results contribute to earthquake hazard assessment, particularly in light of the recent earthquakes in Türkiye, Syria and Morocco. Moreover, by proposing these open-access codes, I advocate for a collaborative effort within the geodetic community to adopt common procedures for cleaning and combining GNSS velocities. This standardisation will enhance the reliability and comparability of geodetic measurements and streamline collaborative efforts and insights across different research groups.

The methodology implemented in the code combines some of the approaches by previous research on aggregated GNSS velocity fields, including [Piña-Valdez., et al., (2022)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021JB023451) and [Zeng et al., (2022)](https://pubs.geoscienceworld.org/ssa/srl/article-abstract/93/6/3121/617675/GPS-Velocity-Field-of-the-Western-United-States).

---

## 2) Overview of the code 

The code is composed of a master Jupyter notebook called **`FIROCO_GNSS.ipynb`**, an input folder called `raw_input` containing the input velocity fields as column-formatted text files with `.raw` extension. A folder called `scripts` contains ancillary Python scripts to filter, combine and rotate GNSS velocity fields. The folder structure is organised as follows:

<pre>
📦FICORO_GNSS  
 ┣ 📂raw_input  
 ┃ ┣ 📜alchalbi_2013.raw  
 ┃ ┣ 📜bahrouni_2020.raw  
 ┃ ┣ 📜billi_2023.raw  
 ┃ ┣ 📜bougrine_2019.raw  
 ┃ ┣ 📜briole_2021.raw  
 ┃ ┣ 📜castro_2021.raw  
 ┃ ┣ 📜england_2016.raw  
 ┃ ┣ 📜ergintav_2023.raw  
 ┃ ┣ 📜euref_all.raw  
 ┃ ┣ 📜euref_ch8.raw  
 ┃ ┣ 📜floyd_2023.raw  
 ┃ ┣ 📜gomez_2020.raw  
 ┃ ┣ 📜graham_2021.raw  
 ┃ ┣ 📜hamiel_2021.raw  
 ┃ ┣ 📜khorrami_2019.raw  
 ┃ ┣ 📜kurt_2023.raw  
 ┃ ┣ 📜mcclusky_2010.raw  
 ┃ ┣ 📜nocquet_2012.raw  
 ┃ ┣ 📜ozdemir_2019.raw  
 ┃ ┣ 📜ozkan_2022.raw  
 ┃ ┣ 📜pinaValdes_2022.raw  
 ┃ ┣ 📜reilinger_2006.raw  
 ┃ ┣ 📜saleh_2015.raw  
 ┃ ┣ 📜serpelloni_2022.raw  
 ┃ ┣ 📜sokhadze_2018.raw  
 ┃ ┣ 📜viltres_2020.raw  
 ┃ ┗ 📜viltres_2022.raw  
 ┣ 📂scripts  
 ┃ ┣ 📜coherence_filter.py  
 ┃ ┣ 📜combine_vel.py  
 ┃ ┣ 📜lognorm_filter.py  
 ┃ ┣ 📜plot_maps_filtering.py  
 ┃ ┗ 📜plot_rotated_vels.py  
 ┣ 📜FICORO_GNSS.ipynb  
 ┗ 📜README.md
</pre>

---

#### FICORO_GNSS code in a nutshell:

**Key Steps:**

1. **Data importing:** Load raw data from an input directory.

2. **Filtering by uncertainty distribution:** Analyse the lognormal distribution of velocity uncertainties, excluding stations that fall outside the 99 percentile in the East and North velocity components.

3. **Z-score Filtering:** Remove stations if velocity magnitudes in the East and North components diverge over 2 sigma from the mean, considering a radius of 20 km.

4. **Velocity field alignment to a common reference frame:** Implement a least squares approach to align all the data sets to a reference velocity field using a 6-parameter Helmert transformation (3 translations and 3 rotations), leveraging on repeated stations in both the input and reference data sets.

5. **Velocity field rotation:** Rotate velocity fields to different reference frames using published Euler vectors.

6. **Velocity field combination:** Combine individual velocity fields into a single comprehensive velocity field.

	1. Station Proximity Analysis: Cluster stations within a 1 km range of one another (they are probably the same station reported by different studies. Although, in some cases, collocated stations might be included). For such proximate stations:
	
		1. Implements the Interquartile Range (IQR) method to detect outliers, omitting solutions displaying disparities in magnitude, azimuthal direction, or both. The thresholds for outlier detection are set as: 	

			- Lower threshold = $Q1 -1.5 * IQR$
			- Upper threshold = $Q3 + 1.5 * IQR$

		2. Calculate median velocities for the East, North, and vertical components. Handle vertical velocities when some of the input solutions do not include verticals.

		3. Derive uncertainties for each velocity component, taking the median value of stations within a $0.01^\circ$ (1.11 km) radius.

		4. Save the combined velocity field data and plot velocity maps.

---

Here is a detailed breakdown of what the master Jupyter Notebook and ancillary scripts do:

### `FICORO_GNSS.IPYNB` 

This script operates as the orchestrator, calling various specialised scripts and organising the flow of data through the filtering, combining, and rotating processes, ultimately aimed at preparing GNSS velocity fields for further geospatial analysis.

This master Jupyter Notebook follows a streamlined process composed of 8 steps:

- Step 1: Define input & output parameters
- Step 2: Load input velocity fields
- Step 3: Identify and remove outliers
- Step 4: Plot filtered GNSS velocities and outliers
- Step 5: Align velocity fields to the International Terrestrial Reference Frame (ITRF14)
- Step 6: Rotate velocity fields using Euler poles
- Step 7: Combine velocity fields using a least squares approach
- Step 8: Plot the combined velocity field in different reference frames

---

#### Step 1: Define input & output parameters

1.  **Library Imports and Path/Parameter Definitions**:
    
    - Import required libraries for handling file systems (`os`), running subprocesses (`subprocess`), working with dates and times (`datetime`), moving files (`shutil`), data processing (`pandas`), scientific computing (`numpy`) and plotting geophysical data (`pygmt`). 
    -   Define paths for input files, output destinations, scripts, and results. These paths are essential for organizing files and directories for the processing that follows. Also, specific parameters related to file types are set.
    
2.  **Input Parameters**:
     
    - Define directory and file structures, including paths for raw inputs, output formatting, scripts, and various stages of processing results.
    
3.  **Script Path Definitions**:
      
    - Paths to ancillary Python scripts are declared. These scripts perform log-normal and coherence filtering, velocity combination, and plotting GNSS velocity maps.

---

#### Step 2: Load input velocity fields

This step handles the initial processing of the input files, ensuring they are available, accessible, and formatted correctly for the subsequent stages of the pipeline. This process is essential for preparing the raw velocity data for further filtering, combination, and analysis.

1.  **Initial status print**:
      
    - The code prints to the console, indicating it is in the process of checking the input files. 
    
2.  **Input path validation**:
      
    - The code first checks whether the input directory exists. If the directory doesn't exist, the script prints an error message and terminates with an exit code of 1, indicating an abnormal termination.
    -   It then checks if the input directory is empty (i.e., no input files). Again, if this check fails, it's considered a fatal error, leading to the script's termination.
    
3.  **Output path handling**:
    
    - The script checks for the existence of the output directory (containing formatted input velocity fields).
    -   If the directory exists, it cleans it by removing all files within. This step ensures that previous run outputs don't mix with the current run.
    -   If the directory doesn't exist, it's created. This setup ensures the script is self-sufficient in managing its required resources.
    
4.  **File loading**:
      
    - The script lists all files in the input directory and begins processing each one that matches the expected file extension (defined in step 1 as `.raw`).
    -   For each valid input file, the code: 
	    * Extracts the base name of the file (name without extension). 
	    * Prints an ongoing status message indicating the translation is being performed. This serves as a log point, helping track each file's processing. 
	    * Reads the content of the file, splitting it into individual lines and further tokenizing each line by spaces, thus forming a list of lists (each representing a row of space-separated values from the file). 
	    * Prepares a header line indicating the 13 columns of the expected format. 
	    * Joins the tokenized lines back into a single string per line and prefixes all the lines with the header, creating a fully formatted text block. 
	    * Writes the formatted content into a new file in the output directory (default: raw_input_column_formatted) with an updated extension (default: `.vel`). This action transforms the file into the '.vel' format expected for subsequent operations.

In summary, the second step is dedicated to validating and preparing the input data, ensuring it's in the right format and the right path for further stages of the pipeline. It avoids involuntary data mixing and provides console logs for monitoring the process.

---

#### Step 3: Identify and remove outliers

This step deals with the cleaning and filtering of the input GNSS velocity fields. The operations use ancillary Python scripts, specifically designed for these tasks. Here are the detailed steps:

1.  **Initial status print**:
      
    - The script provides a console output indicating the initiation of the GNSS velocity fields' cleaning process. 
    
2.  **Lognorm Filtering**:
     
    - The code checks whether the 'lognorm' filter script exists at the specified path. If it does not, an error message is printed, and the program terminates. This check prevents the program from crashing and provides a clear indication of the failure.
    -  If the 'lognorm' filter script exists, the master script prepares to execute it. It defines paths where the input data is stored and where to store the outputs after the 'lognorm' script processes the data. This process includes:
        -   Input data folder (`input_lognorm_path`).
        -   Output folder for filtered data (`output_lognorm_path`).
        -   A directory for data excluded during filtering (`excluded_lognorm_path`).
        -   A folder for figures (`figure_folder_path`).
    -   The code then calls `subprocess.run()`, executing the 'lognorm' filter script with the necessary input and output arguments. 
    
3.  **Coherence Filtering**:
      
    - Similar to the 'lognorm' filter script, the master script checks the existence of the 'coherence' filter script. If the script is missing, the program prints an error message and exits, ensuring that the failure point is known.
    -  If the script is available, it executes the 'coherence' script using `subprocess.run()`, but this time, it passes the output path from the 'lognorm' filtering as an input. The 'coherence' filtering is a subsequent step after 'lognorm' filtering, and it processes the already cleaned data further. 

---

#### Step 4: Plotting GNSS velocities and outliers for each study:

Step 4 in the workflow pertains to the visual representation of the GNSS velocity fields. This step is crucial for several reasons, as it transitions from data cleaning and manipulation to visual interpretation and analysis. Here's what happens during this step:

1.  **Initial status print**:
    
    -   The code begins with printing lines to the console, indicating that the plotting of velocity fields has been initiated. 
    
2.  **Plotting script execution check**:
    
    -   The code checks whether the required plotting script (`plot_maps_filter_path`) exists in the specified path. This is a preventative measure ensuring that the workflow doesn't proceed with a missing file, which would lead to errors.
    -   If the script is found, it is executed using `subprocess.run()`, which is a Python standard library method for invoking subprocesses. In this context, it's used to run a Python script that handles the map plotting using the pygmt module.
    -   If the script is not found, an error message prints, and the program terminates by calling `exit(1)`, indicating that an error has occurred. 
    - For each input velocity field, the plotting script takes the cleaned data and outliers from the previous velocity filtering steps and creates a visual representation on a map showing vectors indicating tectonic plate motions at each GNSS station. Different colours are used to denote whether a station passed the filtering step or was classified as an outlier.

This step allows users to visualise both, cleaned velocities and outliers with different colours on the same map for each of the input velocity fields. Below I show the results for one of the velocity fields by [Billi., et al., (2023)](https://doi.org/10.1016/j.epsl.2022.117906)covering northern Africa. The left panel shows velocities and outliers based on the two criteria discussed in Step 3, the centre and right panels show the the lognormal fit to the velocity uncertainty distribution in East and North velocity components expressed in mm/yr. Stations with velocity uncertainties above the 99th lognormal percentile (black dashed line) in East or North velocity components are classified as outliers. 

![Horizontal and vertical velocity field in the Mediterranean and Middle East areas](Readme_figures/filtering.jpg)

---

#### Step 5: Align velocity fields to the International Terrestrial Reference Frame (ITRF14):

Step 5 is a comprehensive data preparation and processing stage ensuring that GNSS data is correctly formatted for use with the external Fortran script VELROT. This program performs a least squares inversion leveraging common sites between a reference velocity field and an input velocity field, performing a 6-parameter Helmert transformation (3 rotations and 3 translations without scale). The master velocity field is expressed in a no-net-rotation reference frame (ITRF2014), and each input file is aligned with the master velocity field.

Step 5 involves several directories and files, each serving a specific purpose in the data processing workflow. I'll detail each path and its role within the script:

1.  **`coherence_results_path`**:
    
    -  This directory contains the input CSV files that the script processes. These files correspond to the output from the coherence filtering step. The script checks whether the directory is empty before proceeding.
    
2.  **`input_files_rot_path`**:
    
    -  The code uses this directory to store '.vel' files, converted from the original CSVs. If the directory doesn't exist, the script creates it; if it does, the script clears any existing files to prevent data confusion.
    
3.  **`ITRF14_vel_path`**:
    
    -  After converting and processing the '.vel' files, the script uses this directory to handle the outputs that have been aligned with the ITRF14 reference frame. Similar to before, the script prepares this directory by creating or cleaning it as necessary.
    
4.  **`lnk_file_path` and `lnk_folder_path`**:
    
    -  These paths relate to a 'link' file required for the alignment process with the Fortran code VELROT. The script ensures the link file exists, creating it if it doesn't. 
    
5.  **`reference_vel_path`**:
    
    -  This is a critical file acting as a reference for velocity fields. The script checks explicitly that this file exists because it's necessary for the alignment of velocity fields to the ITRF14 standard. In the code, the reference velocity field to which all the other input velocity fields are aligned is that from Serpelloni et al., (2022), which is expressed in ITRF14 and represents one of the most reliable and extensive velocity fields in the data set. A large number of stations in the reference velocity field is crucial because the least squares alignment and estimation of Helmert parameters is done relying on common stations between the input velocity fields and the master velocity field.
    
6.  **`results_path`, `igb_nocomb_path`, and `igb_nocomb_subpath`**:
    
    -  These directories are involved in the final stages of the process, where the script stores the final outputs. Specifically, 'igb14_no_comb' is a structured directory where the script organizes the processed data aligned to the 'igb14' reference frame. The code handles the creation of these directories if they don't exist to properly archive the final products.

---
#### Step 6: Rotate velocity fields using Euler poles

---
#### Step 7: Combine velocity fields

![Number of independent velocity estimates at each GNSS station](Readme_figures/num_estimates.jpg)
---
#### Step 8: Plot combined velocity field in different reference frames

![Horizontal and vertical velocity field in the Mediterranean and Middle East areas](Readme_figures/gps_map.jpg)

---
#### To do: Unit testing

- Check expected number of rows and columns in output files
- Compare combined velocity field with published velocities from the International Terrestrial Reference Frame 2014 (check magnitude of velocity residuals at common GNSS sites). Plot a map showing residual velocities.
