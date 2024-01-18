# Filtering, Combining and Rotating GNSS Velocity Fields in Python (FICORO)

## 1) Introduction

The Alpine-Himalayan belt provides a unique opportunity to investigate plate dynamics and large-scale crustal deformation and earthquake hazard given the complex tectonic interactions involving multiple plates and tectonic regimes. The expansion of regional GNSS networks and the availability of published velocities have facilitated a broader understanding of active tectonics in this region. However, despite these advancements, few attempts have been made to integrate the available GNSS velocities, and consensus on the optimal methods for filtering and harmonizing geodetic datasets remains elusive. The primary objective of this project is to filter, combine and rotate GNSS velocity fields from column-formatted text files. These data are obtained from different geodetic studies of lithospheric deformation, with each input file comprising 13 distinct columns, including 'Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj', 'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig', 'Stat'. The goal is to combine all the published GNSS velocity fields along the Alpine-Himalayan belt, rotate them into different reference frames, filter outliers, and manage repeated stations that might differ in names or have coordinates varying by $\leq 0.01^\circ$ (1.11 km) across different studies. By providing insights into plate motions and crustal deformation, our results contribute to earthquake hazard assessment, particularly in light of the recent earthquakes in Türkiye, Syria and Morocco. Moreover, by proposing these open-access codes, I advocate for a collaborative effort within the geodetic community to adopt common procedures for cleaning and combining GNSS velocities. This standardisation will enhance the reliability and comparability of geodetic measurements and streamline collaborative efforts and insights across different research groups.

The methodology implemented in the code combines some of the approaches by previous research on aggregated GNSS velocity fields, including [Piña-Valdez., et al., (2022)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021JB023451) and [Zeng et al., (2022)](https://pubs.geoscienceworld.org/ssa/srl/article-abstract/93/6/3121/617675/GPS-Velocity-Field-of-the-Western-United-States).

---

**Note:** This code currently relies on the installation of GAMIT/GLOBK, as it utilizes the Fortran codes VELROT and CVFRAME included with GAMIT/GLOBK to align and rotate velocity fields. In upcoming releases, I plan to provide Python scripts that will eliminate the need for GAMIT/GLOBK. Until then, please ensure that GAMIT/GLOBK is installed before running this code.

---

## 2) Overview of the code 

This software comprises a main Jupyter notebook named **`FICORO_GNSS.ipynb`** and an input folder titled `raw_input`. Within the `raw_input` folder, you'll find the input velocity fields stored as column-formatted text files with the `.raw` extension. To facilitate data processing, there's a `scripts` folder containing additional Python scripts designed for filtering, rotating and combining GNSS velocity fields. If you need to manually remove outliers from the data, you can utilize the `manual_filter` folder, which houses a CSV file that enables you to define specific geographic coordinates (latitude and longitude) and corresponding radii (in kilometers) for the removal of outliers from the combined velocity fields in diferent reference frames. To provide clarity, the folder structure is organized as follows:

<pre>
📦FICORO_GNSS
 ┣ 📜FICORO_GNSS.ipynb
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
 ┃ ┣ 📜liang_2013.raw
 ┃ ┣ 📜mcclusky_2010.raw
 ┃ ┣ 📜nocquet_2012.raw
 ┃ ┣ 📜ozdemir_2019.raw
 ┃ ┣ 📜ozkan_2022.raw
 ┃ ┣ 📜pinaValdes_2022.raw
 ┃ ┣ 📜reilinger_2006.raw
 ┃ ┣ 📜saleh_2015.raw
 ┃ ┣ 📜serpelloni_2022.raw
 ┃ ┣ 📜sokhadze_2018.raw
 ┃ ┣ 📜stamps.raw
 ┃ ┣ 📜viltres_2020.raw
 ┃ ┣ 📜viltres_2022.raw
 ┃ ┣ 📜wang_barbot_2023.raw
 ┃ ┣ 📜wang_shen_2020.raw
 ┃ ┗ 📜wedmore_2021.raw
 ┣ 📂scripts
 ┃ ┣ 📜coherence_filter.py
 ┃ ┣ 📜combine_vel.py
 ┃ ┣ 📜lognorm_filter.py
 ┃ ┣ 📜plot_maps_filtering.py
 ┃ ┗ 📜plot_rotated_vels.py
 ┣ 📂manual_filter
 ┃ ┗ 📜filter_criteria.csv
 ┣ 📂Readme_figures
 ┃ ┣ 📜filtering.jpg
 ┃ ┣ 📜gps_map.jpg
 ┃ ┗ 📜num_estimates.jpg
 ┗ 📜README.md
</pre>

---

#### FICORO_GNSS code in a nutshell:

**Key Steps:**

1. **Data importing:** Load raw data from an input directory.

2. **Filtering by uncertainty distribution:** Analyse the lognormal distribution of velocity uncertainties, excluding stations that fall outside the 99 percentile in the East and North velocity uncertainty components.

3. **Z-score Filtering:** Remove stations if velocity magnitudes in the East and North velocity components diverge over 2 sigma from the mean, considering a radius of 20 km. Additionally, the code allows applying less stringent filtering criteria (n sigma) based on the name of the input file, allowing for a customizable approach to data filtering.

4. **Velocity field alignment to a common reference frame:** Implement a least squares approach to align all the data sets to a reference velocity field using a 6-parameter Helmert transformation (3 translations and 3 rotations), leveraging on repeated stations in both the input and reference data sets.

5. **Velocity field rotation:** Rotate velocity fields to different reference frames using published Euler vectors.

6. **Velocity field combination:** Combine individual velocity fields into a single comprehensive velocity field.

	1. Station Proximity Analysis: Cluster stations within a $0.01^\circ$ (1.11 km) range of one another (they are probably the same station reported by different studies. Although, in some cases, collocated stations might be included too). For such proximate stations:
	
		1. Implements the Interquartile Range (IQR) method to detect outliers, omitting solutions displaying disparities in magnitude, azimuthal direction, or both. The thresholds for outlier detection are set as: 	

			- Lower threshold = $Q1 -1.5 * IQR$
			- Upper threshold = $Q3 + 1.5 * IQR$

		2. Calculate median velocities for the East, North, and vertical components. Handle vertical velocities when some of the input solutions do not include verticals.

		3. Derive uncertainties for each velocity component, taking the median value of stations within a $0.01^\circ$ (1.11 km) radius.

		4. Save the combined velocity field data and plot velocity maps.

7. **Manual filtering**: Remove outliers based on geographical coordinates and radii, generating cleaned data files and logs of removed stations. The script uses parallel processing to handle multiple input velocity fields in diferent reference frames efficiently. 

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
- Step 8: Manual filtering of outliers
- Step 9: Plot the combined velocity field in different reference frames

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
    - For each input velocity field, the plotting script takes the cleaned data and outliers from the previous velocity filtering steps and creates a map showing vectors indicating the velocity at each GNSS station. Different colours are used to denote whether a station passed the filtering step or was classified as an outlier.

This step allows users to visualise both, cleaned velocities and outliers with different colours on the same map for each of the input velocity fields. Below I show the results for one of the velocity fields by [Billi., et al., (2023)](https://doi.org/10.1016/j.epsl.2022.117906) covering northern Africa. The left panel shows velocities and outliers based on the two criteria discussed in Step 3, the centre and right panels show the the lognormal fit to the velocity uncertainty distribution in East and North velocity components expressed in mm/yr. Stations with velocity uncertainties above the 99th lognormal percentile (black dashed line) in East or North velocity components are classified as outliers. 

![Horizontal and vertical velocity field in the Mediterranean and Middle East areas](Readme_figures/filtering.jpg)

---

#### Step 5: Align velocity fields to the International Terrestrial Reference Frame (ITRF14):

Step 5 is a comprehensive data preparation and processing stage ensuring that GNSS data is correctly formatted for use with the external Fortran script VELROT, which is distributed with the [GAMIT/GLOBK software package](http://geoweb.mit.edu/gg/). VELROT performs a least squares inversion leveraging common sites between a reference velocity field and an input velocity field, performing a 6-parameter Helmert transformation (3 rotations and 3 translations without scale). The master velocity field is expressed in a no-net-rotation reference frame (ITRF2014), and each input file is aligned with the master velocity field (Serpelloni et al., 2022).

1. **Column Formatting for velrot Compatibility**

- **Folder Checks**:
  - Checks `coherence_results_path` for data; terminates script if empty.
  - Ensures `input_files_rot_path` exists, creates it if necessary.
    - For each `.vel` file in `input_files_rot_path`, the script dynamically creates a new folder in `ITRF14_vel_path`.
    - Each folder is named after the base name of the corresponding velocity file, ensuring a clear link between the input file and its processing directory.

- **Data Formatting and file/folder management**:
  - Processes `.csv` files in `coherence_results_path`, converting them to `.vel` format.
  - Writes formatted data to `input_files_rot_path` for alignment.
  - The script copies three essential types of files into these newly created folders:
  - The input velocity file (`*.vel`).
  - The reference velocity file (`reference_vel`).
  - The `VELROT` link file (`lnk_file`).

This setup guarantees that each folder contains all necessary files for `velrot` alignment, isolated from other datasets.

2. **Aligning Velocity Fields to ITRF14 with Log Creation**

- **Pre-Alignment Setup**:
  - Confirms `input_files_rot_path` is populated.
  - Manages `ITRF14_vel_path` for output data.
  - Creates `lnk_file` if absent.
  - Verifies existence of `reference_vel` file.

- **Conditional Processing to Avoid Redundant Rotation**:
  - Identifies if a `.vel` file is the `reference_vel` (serpelloni et al., 2022 velocity field).
  - Instead of running `velrot`, the script renames the reference velocity file, appending _igb14.vel to its base name.
  - Proceeds with rotation for other files.

- **Rotation and Log File Creation**:
  - Executes the `velrot` command for all the other `.vel` files, aligning them to Serpelloni's velocity field expressed in the IGB14 realisation of ITRF14.
  - Captures output statistics from the alignment process in log files (`*_align.log`).

- **Post-Rotation Processing**:
  - Cleans and extracts data from `velrot` output.
  - Produces aligned `.vel` files with `_igb14.vel` suffix.

- **Organizing Aligned Data**:
  - Copies aligned velocity files into `igb_nocomb_path/igb_nocomb_subpath`.

Step 5 involves several directories and files, each serving a specific purpose in the data processing workflow. I'll detail each path and its role within the script:

-  **`coherence_results_path`**:
    
    -  This directory contains the input CSV files that the script processes. These files correspond to the output from the coherence filtering step. The script checks whether the directory is empty before proceeding.
    
-  **`input_files_rot_path`**:
    
    -  The code uses this directory to store '.vel' files, converted from the original CSVs. If the directory doesn't exist, the script creates it; if it does, the script clears any existing files to prevent data confusion.
    
-  **`ITRF14_vel_path`**:
    
    -  After converting and processing the '.vel' files, the script uses this directory to handle the outputs that have been aligned with the ITRF14 reference frame. Similar to before, the script prepares this directory by creating or cleaning it as necessary.
    
-  **`lnk_file_path` and `lnk_folder_path`**:
    
    -  These paths relate to the VELROT 'link' file required for the alignment process with the Fortran code VELROT, distributed with GAMIT/GLOBK. The script ensures the link file exists, creating it if it doesn't. 
    
-  **`reference_vel_path`**:
    
    -  This is a critical file acting as a reference for velocity fields. The script checks explicitly that this file exists because it's necessary for the alignment of velocity fields to the ITRF14 standard. In the code, the reference velocity field to which all the other input velocity fields are aligned is that from Serpelloni et al., (2022), which is expressed in ITRF14 and represents one of the most reliable and extensive velocity fields in the data set. A large number of stations in the reference velocity field is crucial because the least squares alignment and estimation of Helmert parameters is done relying on common stations between the input velocity fields and the master velocity field.
    
-  **`results_path`, `igb_nocomb_path`, and `igb_nocomb_subpath`**:
    
    -  These directories are involved in the final stages of the process, where the script stores the final outputs. Specifically, 'igb14_no_comb' is a structured directory where the script organizes the processed data aligned to the 'igb14' reference frame. The code handles the creation of these directories if they don't exist to properly archive the final products.

---
#### Step 6: Rotate velocity fields using Euler poles

This step focuses on rotating GNSS velocity fields into different tectonic reference frames using Euler poles.

1. **Setting Up Destination Folders**:
   - Destination folders for each reference frame (`arab`, `eura`, `nubi`, `sina`, `anat`) and `igb14` are set up within `results_path/igb_nocomb_path`.
   
   - **Folder Paths**: The script sets up distinct destination folders for each tectonic plate (`arab`, `eura`, `nubi`, `sina`, `anat`) and for the `igb14` reference frame.
        - **Variables**: 
        - `arab_dest_folder`: Destination folder for velocity fields expressed in Arabia-fixed reference frame
        - `eura_dest_folder`: Destination folder for velocity fields expressed in Eurasia-fixed reference frame
        - `nubi_dest_folder`: Destination folder for velocity fields expressed in Nubia-fixed reference frame
        - `sina_dest_folder`: Destination folder for velocity fields expressed in Southern Sinai-fixed reference frame
        - `anat_dest_folder`: Destination folder for velocity fields expressed in East Anatolia-fixed reference frame
        - `igb14_dest_folder`: Destination folder for for velocity fields expressed in IGB14/ITRF14 reference frame
    - **Folder Creation**: Each destination folder is created using `os.makedirs`, ensuring the script can proceed without interruption due to missing directories.

2. **Defining Euler Poles**:
   - Euler pole parameters for Arabia, Eurasia, Nubia, Sinai, and Anatolia are defined based on recent studies.

3. **Processing Velocity Files**:
   - Iterates through each `.vel` file in `igb14_dest_folder`.
   - Removes `_igb14` suffix from the base name.
   - Applies Euler pole rotation for each tectonic plate using the `cvframe` script distributted with GAMIT/GLOBK.
   - Names rotated files as `"{base}_{ref_frame}.vel"`.

4. **File Management and Rotation**:
   - Moves each rotated file to its respective destination folder based on the reference frame.

5. **Completion Message**:
   - Prints a completion message post-rotation.

---
#### Step 7: Combine velocity fields

Step 7 involves the combination of GNSS velocity fields for each reference frame, using a parallel processing approach.

1. Combining GNSS Velocities
- Function `combine_velocities(ref_frame)`:
  - A function designed to combine GPS velocities for a given reference frame.
  - It utilizes the `combination_script_path` to process data in the specified `ref_frame_folder`.
- **Process**:
  - For each reference frame, the function is called to combine GNSS velocity data in that specific frame.

2. Parallel Execution
- **ThreadPoolExecutor**:
  - The script uses `concurrent.futures.ThreadPoolExecutor` for parallel execution.
  - This approach significantly speeds up the process by handling multiple reference frames simultaneously.
- **Total Workers**:
  - The number of workers is set to the number of reference frames plus one (for IGB14).
- **Execution**:
  - The `executor.map` function is used to apply `combine_velocities` to each reference frame.
  - Additionally, the IGB14/ITRF14 reference frame combination is handled separately.

3. Error Handling
- The script checks for the existence of the `combination_script_path`.
- If the script is not found, an error message is displayed, and the process is terminated.

![Number of independent velocity estimates at each GNSS station](Readme_figures/num_estimates.jpg)

---
#### Step 8: Manual filtering of outliers

In this step, multiple velocity fields in diferent reference frames are filtered in parellel, applying user-defined filtering based on coordinates and radii. The user can edit the filtering criteria file `./manual_filter/filter_criteria.csv`, which is a space-separated CSV with columns center_lon (in degrees), center_lat (in degrees), radius (in kilometers). Stations located within the specified radius around the coordinate provided will be excluded from the final combined velocity field. Clean data and log files listing removed stations are saved in the folder `./results/combined_velocities/manual_filter/`.

---
#### Step 9: Plot combined velocity field in different reference frames

In this step, the code executes the plotting script `plot_rotated_script_path`, corresponds to `scripts/plot_rotated_vels.py`.  
This script plots the horizontal GNSS velocity fields in different reference frames given an input folder containing the velocity fields as CSV files.

1. Function `plot_gps_velocity_fields(folder_path)`:
- **Purpose**: Display maps of GPS velocity fields from CSV files.
- **Process**:
  - **Find CSV Files**: The function searches for all CSV files in the specified input `folder_path`.
  - **Figure Creation**: For each CSV file, a PyGMT figure is created.
    - Sets region, projection, and frame for the map.
    - Adds shaded topography, coastlines, and custom color palettes.
  - **Data Reading**: Reads CSV file data into a pandas DataFrame, setting appropriate column names.
  - **Data Validation**: Skips plotting if the DataFrame is empty.
  - **Velocity Vector Creation**: 
    - Extracts coordinates, velocity components, and uncertainties (which are not used nor plotted in the current version of the code).
    - Calculates the velocity magnitude and normalizes it.
    - Creates a list of vectors for plotting.
  - **Plotting**:
    - Plots GPS velocity vectors with specified style and color.
    - Adds a scale bar to the map.
  - **Displaying Figures**: Shows the plotted figure and prints a status message for each reference frame.

![Horizontal and vertical velocity field in the Mediterranean and Middle East areas](Readme_figures/gps_map.jpg)

---
#### To do: Unit testing

- Check expected number of rows and columns in output files
- Compare combined velocity field with published velocities from the International Terrestrial Reference Frame 2014 (check magnitude of velocity residuals at common GNSS sites). Plot a map showing residual velocities.
