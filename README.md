# FICORO_GNSS

## Python scripts to filter, combine and rotate GNSS velocity fields

This program filters, combines and rotates GNSS velocity fields from column-formatted text files. Each input file should have 13 columns, including 'Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj', 'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig','Stat'. The goal is to combine published velocity fields in different reference frames, filter outliers, and manage repeated stations that might differ in names or have coordinates varying by $\leq 0.01^\circ$ (1.11 km) across different studies.

**Key Steps:**

1. **Data importing:** Load raw data from an input directory.

2. **Filtering by uncertainty distribution:** Analyse the lognorm distribution of velocity uncertainties, excluding stations that fall outside the 99 percentile in the East and North velocity components.

3. **Z-score Filtering:** Remove stations if velocity magnitudes in the East and North component diverge over 2 sigma from the mean, considering a radius of 20 km.

4. **Velocity field alignment to a common reference frame:** Implement a least squares approach to align all the data sets to a reference velocity field using a 6-parameter Helmert transformation (3 translations and 3 rotations), leveraging on repeated stations in both the input and reference data sets.

5. **Velocity field rotation:** Rotate velocity fields to different reference frames using published Euler vectors.

6. **Velocity field combination:** Combine individual velocity fields into a single comprehensive velocity field.

	1. Station Proximity Analysis: Cluster stations within a 1 km range of one another (they are probably the same station reported by different studies. Although, in some cases collocated stations might be included). For such proximate stations:
	
		1. Implements the Interquartile Range (IQR) method to detect outliers, omitting solutions displaying disparities in magnitude, azimuthal direction, or both. The thresholds for outlier detection are set as: 	

			- Lower threshold = $Q1 -1.5 * IQR$
			- Upper threshold = $Q3 + 1.5 * IQR$

		2. Calculate median velocities for the East, North, and vertical components. Handle vertical velocities when some of the input solutions do not include verticals.

		3. Derive uncertainties for each velocity component, taking the median value of stations within a $0.01^\circ$ (1.11 km) radius.

		4. Save the combined velocity field data and plot velocity maps.
