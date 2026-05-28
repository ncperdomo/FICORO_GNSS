
# FICORO_GNSS

[![Language](https://img.shields.io/badge/python-3%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](https://github.com/ncperdomo/FICORO_GNSS/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/682602490.svg)](https://doi.org/10.5281/zenodo.13921188)
[![nbviewer](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.org/github/ncperdomo/FICORO_GNSS/blob/main/FICORO_GNSS.ipynb)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ncperdomo/FICORO_GNSS_binder.git/main?urlpath=%2Fdoc%2Ftree%2FFICORO_GNSS.ipynb)

## An open-source Python software package for filtering, combining and rotating GNSS velocity fields

1. [Introduction](#1-introduction)
2. [Overview of the code](#2-overview-of-the-code)
3. [Usage](#3-usage)
4. [Example outputs](#4-example-outputs)
5. [How to cite](#5-how-to-cite)
6. [License](#6-license)

## 1) Introduction

The expansion of regional GNSS networks and the availability of published velocities have significantly enhanced our understanding of active tectonics. However, despite these advancements, few attempts have been made to integrate the available GNSS velocities at continental and global scales (e.g., [Nocquet J.M, 2012](https://doi.org/10.1016/j.tecto.2012.03.037); [Kreemer et al., 2014]( https://doi.org/10.1002/2014GC005407); [Graham et al., 2018](https://doi.org/10.1029/2017GC007391); [Piña-Valdez., et al., (2022)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021JB023451); [Zeng et al., (2022)](https://pubs.geoscienceworld.org/ssa/srl/article-abstract/93/6/3121/617675/GPS-Velocity-Field-of-the-Western-United-States)), and integrating multiple GNSS velocity fields remains challenging due to the lack of standardised methods for filtering and harmonising these data sets.

FICORO_GNSS is a Python-based suite of tools designed to filter, combine, and rotate multiple GNSS velocity fields into a comprehensive and self-consistent dataset. Each input GNSS velocity file adheres to the GAMIT\GLOBK velocity file format, containing 13 columns: 
```
Lon Lat E.vel N.vel E.adj N.adj E.sig N.sig Corr U.vel U.adj U.sig Stat 
```

The methodology implemented in the code combines some of the approaches by previous research on combined GNSS velocity fields, including [Piña-Valdez., et al., (2022)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021JB023451) and [Zeng et al., (2022)](https://pubs.geoscienceworld.org/ssa/srl/article-abstract/93/6/3121/617675/GPS-Velocity-Field-of-the-Western-United-States).

---

**Note:** Earlier versions of this software (v1.x) depended on the installation of GAMIT/GLOBK, as they used the Fortran programs `VELROT` and `CVFRAME` included in the [GAMIT/GLOBK package](https://geoweb.mit.edu/gg/) for velocity field alignment and rotation. The current release (v2.x) replaces these dependencies with fully integrated Python implementations ([PyVelrot](https://github.com/ncperdomo/PyVelrot) and [PyCvframe](https://github.com/ncperdomo/PyCvframe)), eliminating the need to install GAMIT/GLOBK while preserving the original alignment and rotation functionality. 

## 2) Overview of the code

This software comprises a main Jupyter notebook named **`FICORO_GNSS.ipynb`** and an input folder titled `raw_input`. Within the `raw_input` folder, you'll find the input velocity fields stored as column-formatted text files with the `.raw` extension. To facilitate data processing, there's a `scripts` folder containing additional Python scripts designed for filtering and combining GNSS velocity fields. If you need to manually remove outliers from the data, you can utilise the `manual_filter` folder, which houses a CSV file that enables you to define specific geographic coordinates (latitude and longitude) and corresponding radii (in kilometers) for the removal of outliers from the combined velocity field. 

The folder structure is organised as follows:

<pre>
📦FICORO_GNSS
 ┣ 📜FICORO_GNSS.ipynb
 ┣ 📜README.md
 ┣ 📂scripts
 ┃ ┣ 📜coherence_filter.py
 ┃ ┣ 📜combine_vel.py
 ┃ ┣ 📜frame_registry.py 
 ┃ ┣ 📜lognorm_filter.py
 ┃ ┣ 📜plot_maps_filtering.py
 ┃ ┣ 📜plot_rotated_vels.py
 ┃ ┣ 📜pycvframe.py.py
 ┃ ┣ 📜pyvelrot.py
 ┃ ┗ 📜uncertainty_filter_verticals.py
 ┃ ┗ 📜uncertainty_scaling_combined.py
 ┣ 📂manual_filter
 ┃ ┗ 📜filter_criteria.csv
 ┣ 📂raw_input
 ┃ ┣ 📂levelling_verticals
 ┃ ┃ ┗ 📜wu_2022.raw
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
 ┃ ┣ 📜jolivet_2023.raw
 ┃ ┣ 📜kadirov_rot_karakhanyan_2013.raw
 ┃ ┣ 📜karakhanyan_2013.raw
 ┃ ┣ 📜khorrami_2019.raw
 ┃ ┣ 📜kurt_2023.raw
 ┃ ┣ 📜li_2024.raw
 ┃ ┣ 📜liang_2013.raw
 ┃ ┣ 📜mcclusky_2010.raw
 ┃ ┣ 📜nocquet_2012.raw
 ┃ ┣ 📜ozarpaci_2020.raw
 ┃ ┣ 📜ozbey_2024.raw
 ┃ ┣ 📜ozdemir_2019.raw
 ┃ ┣ 📜ozkan_2022.raw
 ┃ ┣ 📜perry_2018.raw
 ┃ ┣ 📜pinaValdes_2022.raw
 ┃ ┣ 📜reilinger_2006.raw
 ┃ ┣ 📜saleh_2015.raw
 ┃ ┣ 📜serpelloni_2022.raw
 ┃ ┣ 📜sokhadze_2018.raw
 ┃ ┣ 📜stamps_2018.raw
 ┃ ┣ 📜tatar_2012.raw
 ┃ ┣ 📜viltres_2020.raw
 ┃ ┣ 📜viltres_2022.raw
 ┃ ┣ 📜wang_barbot_2023.raw
 ┃ ┣ 📜wang_shen_2020.raw
 ┃ ┣ 📜wedmore_2021.raw
 ┃ ┣ 📜zheng_2017.raw
 ┃ ┗ 📜zubovich_2010.raw
</pre>

---

### FICORO_GNSS in a nutshell:

**Key Steps:**

1. **Filtering stations affected by postseismic transient motions:** Remove stations identified as being affected by postseismic transient motions

2. **Filtering by uncertainty distribution:** GNSS stations with velocity uncertainties exceeding the 99th percentile of the scaled log-normal distribution are removed from input velocity fields, following the approach by [Piña-Valdez., et al., (2022)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021JB023451).

3. **Filtering based on spatial coherence of velocity magnitudes:** Remove stations if velocity magnitudes in the East and North velocity components diverge over 2 sigma from the mean, considering a radius of 20 km. Additionally, the code allows applying geographic-based stringency levels (n-sigma), allowing for a customisable approach to data filtering.

4. **Velocity field alignment to a common reference frame:** Implement a least squares approach to align all the data sets to a reference velocity field using a 6-parameter Helmert transformation (3 translations and 3 rotations), leveraging on repeated stations in both the input and reference data sets.

5. **Velocity field rotation:** Rotate velocity fields to different reference frames using published Euler poles.

6. **Filtering by velocity magnitude and azimuthal direction:** Implement the Interquartile Range (IQR) method to detect outliers, omitting solutions displaying disparities in magnitude, azimuthal direction, or both. The thresholds for outlier detection are set as: 	

	- Lower threshold = $Q1 -1.5 * IQR$
	- Upper threshold = $Q3 + 1.5 * IQR$

7. **Velocity field combination:** Estimate median velocities and uncertainties for the East and North velocity components at collocated stations.

8. **Manual filtering**: Remove outliers based on geographical coordinates and radii, generating cleaned data sets and logs of removed stations. 

9. **Scaling velocity uncertainties:** Horizontal velocity uncertainties are adjusted to match the same percentile in a subjectively chosen target log-normal distribution, following the approach by [Piña-Valdez., et al., (2022)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021JB023451).

10. **Final filtering based on uncertainty distribution:** Stations with velocity uncertainties exceeding the 99th percentile of the scaled log-normal distribution are removed from the final combined velocity field.

## 3) Usage

### Prerequisites

- **Python:** Version 3.7 or higher
- **Python Libraries:** numpy, scipy, matplotlib, pygmt, jupyter, pandas, os, subprocess, datetime, sys, glob, json, time, concurrent, argparse, itertools

### Steps

1. **Clone the Repository:**

    ```bash
    git clone https://github.com/ncperdomo/FICORO_GNSS.git
    cd FICORO_GNSS
    ```

2. **Create a virtual environment (optional but recommended):**
    
    ```bash
    python3 -m venv ficoro_venv
    source ficoro_venv/bin/activate 
    ```

3. **Launch Jupyter Notebook:**
    ```bash
    jupyter notebook
    ```

4. **Open and run FICORO_GNSS.ipynb:**

- Navigate to the repository directory in the Jupyter interface and open the main notebook to start processing GNSS velocity fields.

- Input Data: Place your .raw GNSS velocity files in the `raw_input/` folder. Each file should follow the GAMIT velocity file format. E-N-U adjustments and E-N velocity correlation values can be set to zero if not available.

- Manual Filtering: Use the `manual_filter/` folder to define specific geographic coordinates and radii for outlier removal. Modify the provided CSV file to specify the criteria.


## 4) Example outputs:

**FICORO_GNSS** includes a set of input velocity fields for the Alpine-Himalayan region. By running the main Jupyter notebook, you can filter, rotate, and combine these datasets to generate publication-quality maps using [PyGMT](https://www.pygmt.org/). Below are examples of the outputs you can expect:

![Horizontal and vertical velocity field in the Mediterranean and Middle East areas](Readme_figures/gps_map.png)

## 5) How to Cite

I am committed to promoting reproducibility and open access in geodetic and tectonic research. Please cite FICORO_GNSS according to how the software or associated data products are used in your work.

### A) If you use the FICORO_GNSS software

Please cite the Zenodo software release if you:

- use the filtering, combination, alignment, or rotation workflows implemented in FICORO_GNSS,
- adapt or modify the code,
- incorporate parts of the methodology into another workflow, or
- use the software for scientific analyses or derived products.

**Citation:**

```text
Castro-Perdomo, N. (2024). FICORO_GNSS: An open-source Python software package for filtering, combining and rotating GNSS velocity fields (Version v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.13921189
```

**BibTeX:**

```bibtex
@software{castro_perdomo_ficoro_gnss_2024,
  author       = {Nicolás Castro-Perdomo},
  title        = {{FICORO\_GNSS: An open-source Python software package for filtering, combining and rotating GNSS velocity fields}},
  year         = 2024,
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.13921189},
  url          = {https://doi.org/10.5281/zenodo.13921189}
}
```

---

### B) If you use the combined Alpine–Himalayan GNSS velocity field distributed with FICORO_GNSS

Please cite the accompanying *JGR: Solid Earth* [article](https://doi.org/10.1029/2025JB031738) if you use:

- the combined GNSS velocity fields provided with this repository,
- derived strain-rate products,
- regional subsets of the compiled velocity database, or
- analyses based on the Alpine–Himalayan GNSS compilation.

**Citation:**

```text
Castro-Perdomo, N., Jónsson, S., Klinger, Y., Masson, F., Becker, T. W., & Johnson, K. (2025). Strain rates along the Alpine-Himalayan Belt from a comprehensive GNSS velocity field. Journal of Geophysical Research: Solid Earth, 130(12), e2025JB031738. https://doi.org/10.1029/2025JB031738
```

**BibTeX:**

```bibtex
@article{castro_perdomo_2025_JGR,
  author = {Castro-Perdomo, N. and Jónsson, S. and Klinger, Y. and Masson, F. and Becker, T. W. and Johnson, K.},
  title = {Strain Rates Along the Alpine-Himalayan Belt From a Comprehensive GNSS Velocity Field},
  journal = {Journal of Geophysical Research: Solid Earth},
  volume = {130},
  number = {12},
  pages = {e2025JB031738},
  keywords = {strain rates, active tectonics, crustal deformation, Alpine-Himalayan belt, GNSS, velocity field},
  doi = {https://doi.org/10.1029/2025JB031738},
  url = {https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2025JB031738},
  eprint = {https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2025JB031738},
  note = {e2025JB031738 2025JB031738},
  year = {2025}
}
```

## 6) License

This project is licensed under the MIT License.

---
