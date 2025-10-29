
READ ME File For "[Dataset] Detection of Magnetic Reconnection between A Transpolar Arc and the Interplanetary Magnetic Field"

Dataset DOI: 10.5258/SOTON/D3194 (or 10.5281/zenodo.17077827)

Date that the file was created: October, 2025

-------------------
GENERAL INFORMATION
-------------------

ReadMe Author: NAWAPAT KAWEEYANUN, University of Southampton, ORCID ID: 0000-0003-0634-0164

Date of data collection: 7 August 2024 - 24 October 2025

Information about geographic location of data collection: Southampton, U.K.

Related projects:
None

--------------------------
SHARING/ACCESS INFORMATION
-------------------------- 

Licenses/restrictions placed on the data, or limitations of reuse:

Recommended citation for the data:

This dataset supports the publication:
AUTHORS: Nawapat Kaweeyanun, Robert C. Fear, Betty S. Lanchester, Daniel K. Whiter, Imogen L. Gingell, Andrew N. Fazakerley, Iannis Dandouras, Stephen B. Mende, Larry J. Paxton  
TITLE: Detection of Magnetic Reconnection between A Transpolar Arc and the Interplanetary Magnetic Field
JOURNAL: Geophysical Research Letters
PAPER DOI IF KNOWN:

Links to other publicly accessible locations of the data: N/A

Links/relationships to ancillary or related data sets: N/A


--------------------
DATA & FILE OVERVIEW
--------------------

This dataset contains:

README.txt: this guideline document

Analysis_Codes.zip: compressed folder containing Python files used to conduct analysis

Subfolder files:
- Tsyganenko (Larquier, 2012; Coxon, 2020): directory containing Tsyganenko-96 (T96) magnetic field tracing module, with Python wrapper by Sebastian de Larquier (2012) and John Coxon (2020)
	Subfolder files:
	- .github: directory containing GitHub-related request files (bug_report.md and issue_request.md under ISSUE_TEMPLATE directory and pull_request_template.md)
	- notebooks: directory containing movie illustration of Dungey Cycle (Dungey Cycle.mov) and Python notebooks for T96 model (Introductory Slides.ipynb, Tsyganenko - Python examples.ipynb)
	- original-fortran: directory containing original Fortran code of the T96 model (Geopack-2008_2010-11-30.for, Geopack-2008_2020-01-01.for, Makefile, T02.f, T96.f, test_geopack.f90, test_geopack_py2f.f90)
	- tests: directory containing script to test T96 model (test.py)
	- tsyganenko: directory containing T96 modern operation module in Fortran and Python (__init__.py, convert.py, geopack_py.for, geopack_py.o, geopack_py.pyf, Makefile, T02.f, T02.o, T89.f, T96.f, T96.o, T96_geopack_08.f, trace.py, TS05.f)  	
	- .gitignore: file that specifies intentionally untracked files to ignore
	- LICENSE.md: Copyright statement
	- README.md: Instructions to install tsyganenko operation module within local Python framework
	- setup.py: A script to install tsyganenko operation module within local Python framework
- CDAWeb_Extract.py: Python module to extract and plot data from NASA's Coordinated Data Analysis Web (including OMNI)
- Cluster_Position_Comp_Bplane: Python script to compare positions of two Cluster spacecraft on a plane defined by magnetic field measurement of the first.
- CSA_Extract.py: Python module to extract data from ESA's Cluster Science Archive
- Cluster_Instr_Data.py: Python module to read CDF Cluster data files from Cluster Science Archive
- CSA Organiser.py: Python function to shift downloaded CDF data files from Cluster Science Archive
- Footprint_Aurora_Overlay.py: Python script to overlay Cluster's magnetic footprint on IMAGE's aurora figure.
- FluxLineSum.py: Python module to plot combined Cluster particle fluxes above a defined energy threshold.
- GeotailExtract.py: Python module to extract data from JAXA/DARTS's Geotail webpage.
- IMAGE_Data.py: Python module to read CDF IMAGE data files from Cluster Science Archive.
- IMAGE_Timeseries.py: Python function to generate a timeseries plot of IMAGE's aurora observations.
- IMSG_Overlay.py: Python function to overlay SSUSI/GUVI scan on a grayscale IMAGE output.
- OverlayPlotter.py: Python function to plot Cluster's magnetic footprint on IMAGE's aurora figure.
- PanelPlotter.py: Python module to plot Cluster data as panel plots between specified datetimes.
- SearchPlotter.py: Python function to identify correct IMAGE data file to plot.
- SpeciesComposition.py: Python module to plot ratio between O+ and H+ ion differential energy fluxes.
- SSUSI_GUVI_Data: Python module to extract and plot SSUSI and GUVI scans.
- T96Analyze.py: Python module to extract T96-traced magnetic footprint data and convert its coordinates.
- T96Script.py: Python script to trace Cluster's magnetic footprint using the T96 model (tsyganenko module)
- WalenTest.py: Python module to conduct Walen test between specified datetimes.


Fig2: directory containing spreadsheet files for variables in Fig. 2 of manuscript.

Subfolder files
- C1_Footprint_2002-03-18.csv: Cluster-1's magnetic footprint for 2002-03-18 at 14:30, 14:44, 14:54, and 15:09 UT
- Fig2{Letter}_IMAGE_HighRes_Data_{Instrument}_{Datetime}_{variable}.csv
	Letter = label corresponding to the subplot (A to L)
	Instrument = IMAGE instrument (WIC, S13 or S12)
	Datetime = datetime of data (14:30, 14:44, 14:54, and 15:09 UT on 18 March 2002)
	Variable = whether data is intensity ("image"), magnetic latitude ("mlat") or magnetic local time ("mlt")
- Fig2M_GUVI_Orb1494_Channel3+4_2002077143307_2002077161022_IntN.csv: GUVI radiance data for Orbit 1494 (northern polar cap, used in Fig. 2M).
- Fig2M_GUVI_Orb1494_Channel3+4_2002077143307_2002077161022_MAGLat: GUVI AACGM latitude data for Orbit 1494 (northern polar cap, used in Fig. 2M).
- Fig2M_GUVI_Orb1494_Channel3+4_2002077143307_2002077161022_MLT: GUVI AACGM MLT data for Orbit 1494 (northern polar cap, used in Fig. 2M).
- TIMED_Geo_18Mar02.csv: TIMED trajectory in geographic coordinates between 14:45-15:15 UT on 18 March 2002 (used in Fig. 2M)


Fig3: directory containing spreadsheet files for variables in Fig. 3 of manuscript.

Subfolder files:
- Fig3A_Bvec_OMNI_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv: 
	Magnetic field data from OMNI dataset between 14:15-15:15 UT, 18 March 2002 (used in Fig. 3A).
- Fig3B_ion_fluxU_C1_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv: 
	Ion differential energy flux, sorted by energy and datetime, from Cluster-1's ion spectrometer between 14:15-15:15 UT, 18 March 2002 (used in Fig. 3B)
- Fig3C_e_fluxU_C1_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv: 
	Electron differential energy flux, sorted by energy and datetime, from Cluster-1's electron spectrometer between 14:15-15:15 UT, 18 March 2002 (used in Fig. 3C)
- Fig3C_SC_Pot_C1_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Spacecraft potential from Cluster-1's electric field and wave instrument between 14:15-15:15 UT, 18 March 2002 (used in Fig. 3C)
- Fig3D_Bvec_C1_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv: 
	Magnetic field data from Cluster-1's fluxgate magnetometer between 14:15-15:15 UT, 18 March 2002 (used in Fig. 3D).
- Fig3E_ion_vel_C1_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Ion velocity data from Cluster-1's ion spectrometer between 14:15-15:15 UT, 18 March 2002 (used in Fig. 3E).
- Fig3F_ion_fluxU_GEOTAIL_2002-03-18T14-15-12Z_2002-03-18T15-14-55Z.csv:
	Ion differential flux, sorted by energy and datetime, from Geotail between 14:15-15:15 UT, 18 March 2002 (used in Fig. 3F).



Fig4: directory containing spreadsheet files for variables in Fig. 4 of manuscript.

Subfolder files:
- Fig4A_H+_fluxU_C4_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Hdyrogen ion differential energy flux, sorted by energy and datetime, from Cluster-4's ion spectrometer between 14:15-15:15 UT, 18 March 2002 (used in Fig. 4A)
- Fig4B_He+_fluxU_C4_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Helium ion differential energy flux, sorted by energy and datetime, from Cluster-4's ion spectrometer between 14:15-15:15 UT, 18 March 2002 (used in Fig. 4B)
- Fig4C_O+_fluxU_C4_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Oxygen ion differential energy flux, sorted by energy and datetime, from Cluster-4's ion spectrometer between 14:15-15:15 UT, 18 March 2002 (used in Fig. 4C)
- Fig4D_ObyH_C4_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Ratio between oxygen and hydrogen ion fluxes, sorted by energy and datetime, from Cluster-4's ion spectrometer between 14:15-15:15 UT, 18 March 2002 (used in Fig. 4D)
- Fig4E_H+_fluxPA_C4_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Hydrogen ion differential energy flux, sorted by pitch angle and datetime, from Cluster-4's ion spectrometer between 14:15-15:15 UT, 18 March 2002 (used in Fig. 4E)
- Fig4F_He+_fluxPA_C4_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Helium ion differential energy flux, sorted by pitch angle and datetime, from Cluster-4's ion spectrometer between 14:15-15:15 UT, 18 March 2002 (used in Fig. 4F)
- Fig4G_O+_fluxPA_C4_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Oxygen ion differential energy flux, sorted by pitch angle and datetime, from Cluster-4's ion spectrometer between 14:15-15:15 UT, 18 March 2002 (used in Fig. 4G)
- Fig4H_e_fluxPA_C4_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Electron differential energy flux, sorted by pitch angle and datetime, from Cluster-4's electron spectrometer between 14:15-15:15 UT, 18 March 2002 (used in Fig. 4H)
- Fig4I_e_fluxU180_C4_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Electron differential energy flux, sorted by energy and datetime, at 180 deg pitch angle (field-antiparallel) from Cluster-4's electron spectrometer between 14:15-15:15 UT, 18 March 2002 (used in Fig. 4I)
- Fig4I_SC_Pot_C4_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Spacecraft potential from Cluster-4's electric field and wave instrument between 14:15-15:15 UT, 18 March 2002 (used in Fig. 4I)


FigS1: directory containing spreadsheet files for variables in Fig. S1 of manuscript.

Subfolder files:
- Walen_C1_2002-03-18T14-57-25Z_2002-03-18T15-03-06Z.csv: Data used for Walen test plot (Fig. S1)


FigS2: directory containing spreadsheet files for variables in Fig. S2 of manuscript.

Subfolder files:
- FigS2A_e_fluxU_C1_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Electron differential energy flux, sorted by energy and datetime, from Cluster-1's electron spectrometer between 14:15-15:15 UT, 18 March 2002 (used in Fig. S2A, same as Fig. 3C)
- FigS2A_SC_Pot_C1_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Spacecraft potential from Cluster-1's electric field and wave instrument between 14:15-15:15 UT, 18 March 2002 (used in Fig. S2A, same as Fig. 3C)
- FigS2B_e_fluxU_C2_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Electron differential energy flux, sorted by energy and datetime, from Cluster-2's electron spectrometer between 14:15-15:15 UT, 18 March 2002 (used in Fig. S2B)
- FigS2B_SC_Pot_C2_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Spacecraft potential from Cluster-2's electric field and wave instrument between 14:15-15:15 UT, 18 March 2002 (used in Fig. S2B)
- FigS2C_e_fluxU_C3_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Electron differential energy flux, sorted by energy and datetime, from Cluster-3's electron spectrometer between 14:15-15:15 UT, 18 March 2002 (used in Fig. S2C)
- FigS2C_SC_Pot_C3_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Spacecraft potential from Cluster-3's electric field and wave instrument between 14:15-15:15 UT, 18 March 2002 (used in Fig. S2C)
- FigS2D_e_fluxU_C4_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Electron differential energy flux, sorted by energy and datetime, from Cluster-4's electron spectrometer between 14:15-15:15 UT, 18 March 2002 (used in Fig. S2D)
- FigS2D_SC_Pot_C4_2002-03-18T14-15-00Z_2002-03-18T15-15-00Z.csv:
	Spacecraft potential from Cluster-4's electric field and wave instrument between 14:15-15:15 UT, 18 March 2002 (used in Fig. S2D, same as Fig. 4I)
 

FigS3: directory containing spreadsheet files for variables in Fig. S3 of manuscript.

Subfolder files:
- FigS3A_H+_fluxU_Sun_C4_2002-03-18T14-49-00Z_2002-03-18T14-54-00Z.csv:
	Proton differential energy flux travelling in sunward direction, sorted by energy and datetime, from Cluster-4's ion spectrometer between 14:49 - 14:54 UT, 18 March 2002 (used in Fig. S3A)
- FigS3B_H+_fluxU_Dawn_C4_2002-03-18T14-49-00Z_2002-03-18T14-54-00Z.csv:
	Proton differential energy flux travelling in dawnward direction, sorted by energy and datetime, from Cluster-4's ion spectrometer between 14:49 - 14:54 UT, 18 March 2002 (used in Fig. S3B)
- FigS3C_H+_fluxU_Dusk_C4_2002-03-18T14-49-00Z_2002-03-18T14-54-00Z.csv:
	Proton differential energy flux travelling in duskward direction, sorted by energy and datetime, from Cluster-4's ion spectrometer between 14:49 - 14:54 UT, 18 March 2002 (used in Fig. S3C)
- FigS3D_H+_fluxU_Antisun_C4_2002-03-18T14-49-00Z_2002-03-18T14-54-00Z.csv: 
	Proton differential energy flux travelling in antisunward direction, sorted by energy and datetime, from Cluster-4's ion spectrometer between 14:49 - 14:54 UT, 18 March 2002 (used in Fig. S3D)
- FigS3E_LineSum_C4_2002-03-18T14-49-00Z_2002-03-18T14-54-00Z.csv:
	Sums of proton differential energy flux above 10 keV in sunward, dawnward, duskward, and antisunward direction, from Cluster-4's ion spectrometer between 14:49 - 14:54 UT, 18 March 2002 (used in Fig. S3E)


Relationship between files, if important for context: Analysis code files form a singular Python module that must be placed in the same folder when used for analysis.  

Additional related data collected that was not included in the current data package: N/A

If data was derived from another source, list source: N/A

If there are there multiple versions of the dataset, list the file updated, when and why update was made: N/A


--------------------------
METHODOLOGICAL INFORMATION
--------------------------

Description of methods used for collection/generation of data: 

Cluster and IMAGE data are obtained from ESA's Cluster Science Archive (https://csa.esac.esa.int/). GUVI data are downloaded from the official webpage (https://guvitimed.jhuapl.edu), though this portal will become defunct on 30 September 2025. OMNI data are taken from NASA's Coordinated Data Analysis Web (https://cdaweb.gsfc.nasa.gov/index.html) and JAXA's DARTS portal (https://darts.isas.jaxa.jp/stp/geotail/) respectively.

Script to download Cluster/IMAGE data is available in CSA_Extract.py, though non-footprint Cluster data can be simultaneously obtained and plotted in PanelPlotter.py. Both OMNI and Geotail data can be obtained and plotted in CDAWeb_Extract.py and Geotail_Extract.py. GUVI data files must be downloaded manually from the linked portal above.


Methods for processing the data: 

Raw data is entered into the code files in the Analysis_Code.zip. Follow below instructions to generate Fig. 2-4 and Fig. S1-S3 (information also available in the script). Note that PanelPlotter.py save into a "ClusterPanels" directory while other scripts save into the current directory.

- Fig 2:
    1. Obtain Cluster's magnetometer data via CSA_Extract.py.
    2. Rewrite Cluster's magnetometer data from CDF to CSV via Cluster_Instr_Data.py.
    3. Ensure that directory containing all CSV rewritten files is in the same directory as T96Script.py.
    4. Set up and run T96Script.py to obtain Cluster's magnetic footprint.
    5. Ensure that directory containing all magnetic footprint files (CSV) is in the same directory as Footprint_Aurora_Overlay.py.
    6. (Fig. 2A-2H) Set up and run Footprint_Aurora_Overlay.py to obtain IMAGE WIC and S12 data.
    7. (Fig. 2I) Ensure that relevant GUVI data files are stored in a directory with structure defined in IMSG_Overlay.py.
    8. Set up and run IMSG_Overlay.py.

    Note: T96Script.p, Footprint_Aurora_Overlay.py and IMSG_Overlay.py are already configured to produce Fig. 2, but datetimes and intervals can be adjusted for other explorations.

- Fig 3:
    1. (Fig. 3A) Run script at the end of CDAWeb_Extract.py to obtain OMNI plots.
    2. In PanelPlotter.py, select the parameter list ("paralist") marked Fig. 3 and set the spacecraft list ("sclist") to 'C1' in the main calling script. Set sc_Comp option to False.
    3. (Fig. 3B-3E) Run PanelPlotter.py to obtain Cluster plots.
    4. (Fig. 3F) Ensure that the Geotail data file (not directory) is in the same directory as Geotail_Extract.py.
    5. Run Geotail_Extract.py.

    Note: Datetimes and parameters can be adjusted for other explorations.

- Fig. 4:
    1. (Fig. 4A-4C) In PanelPlotter.py, select the parameter list ("paralist") marked Fig. 4A-4C and set the spacecraft list ("sclist") to 'C4' in the main calling script, and run the script.
    2. (Fig. 4D) Run SpeciesComposition.py with "ObyH" in the parameter list.
    3. (Fig. 4E-4I) In PanelPlotter.py, select the parameter list marked "Fig. 4E-4I" and run the script.

    Note: Datetimes and parameters can be adjusted for other explorations.

- Fig. S1:
    1. Run script at the end of WalenTest.py at pre-defined time intervals to obtain Fig. S1.
    
    Note: Datetimes and parameters can be adjusted for other explorations.

- Fig. S2:
    1. In PanelPlotter.py, select the parameter list ('paralist') marked Fig. S2 and set the spacecraft list ('sclist') to ['C1','C2','C3','C4'] in the main calling script. Set sc_comp option to True.
    2. Run PanelPlotter.py to obtain Fig. S2A-S2D.

- Fig. S3:
    1. Run script at the end of FluxLineSum.py at pre-defined time intervals to obtain Fig. S1.


Software- or Instrument-specific information needed to interpret the data, including software and hardware version numbers: 

- Python 3.7 or above.
- Selected scientific Python modules (e.g., cdflib, cdasws, aacgmv2, check prerequisites listed at start of Python files)
- Key Graphical Products tool, Cluster Science Archive (https://csa.esac.esa.int/csa-web/#graph) for Fig. S2.

Standards and calibration information, if appropriate: N/A

Environmental/experimental conditions: N/A

Describe any quality-assurance procedures performed on the data:

- Code is checked against known examples. Fear et al., (2014) and Fryer et al., (2021) for Cluster data. Phan et al., (2003) for IMAGE data and Walen Test.

People involved with sample collection, processing, analysis and/or submission:

Nawapat Kaweeyanun, Robert C. Fear, Betty S. Lanchester, Daniel K. Whiter, Imogen L. Gingell, Andrew N. Fazakerley, Iannis Dandouras, Stephen B. Mende

--------------------------
DATA-SPECIFIC INFORMATION <Create sections for each datafile or set, as appropriate>
--------------------------

For C1 Footprint data file
Number of variables: 5
Number of cases/rows: 5 columns, 4 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used: datetime, latitude, magnetic local time
1st column: Datetime of footprint (in DD-MM-YYYY hh:mm:ss format)
2nd column: C1 footprint in the ith component of Fig. 2 (km).
3rd column: C1 footprint in the jth component of Fig. 2 (km).
4th column: C1 footprint in the AACGM coordinate latitude (deg).
5th column: C1 footprint in the AACGM coordinate magnetic local time (dimensionless).
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A


For each IMAGE WIC/S13/S12 intensity ("image") data file
Number of variables: 1
Number of cases/rows: 256 columns, 256 rows (WIC); 128 rows, 128 columns (S13/S12) 
All: IMAGE intensity data (R) corresponding to magnetic latitude and local time (in separate files)

For each IMAGE WIC/S13/S12 magnetic latitude ("mlat") data file
Number of variables: 1
Number of cases/rows: 256 columns, 256 rows (WIC); 128 rows, 128 columns (S13/S12) 
All: IMAGE magnetic latitude (degree) corresponding to magnetic local time and intensity (in separate files)

For each IMAGE WIC/S13/S12 magnetic local time ("mlt") data file
Number of variables: 1
Number of cases/rows: 256 columns, 256 rows (WIC); 128 rows, 128 columns (S13/S12) 
All: IMAGE magnetic local time (0-24) corresponding to magnetic latitude and intensity (in separate files)

For GUVI IntN Data file
Number of variables: 1
Number of cases/rows: 318 columns, 318 rows
Variable list:
All = Auroral radiance in the northern polar cap observed along (vertical) and across (horizontal) GUVI's trajectory.

For GUVI MagLat Data file
Number of variables: 1
Number of cases/rows: 318 columns, 318 rows
Variable list:
All = AACGM coordinate latitude at which the radiance is observed in the northern polar cap along (vertical) and across (horizontal) GUVI's trajectory.

For GUVI MLT Data file
Number of variables: 1
Number of cases/rows: 318 columns, 318 rows
Variable list:
All = AACGM coordinate magnetic local time at which the radiance is observed in the northern polar cap along (vertical) and across (horizontal) GUVI's trajectory.

Fig 3A (Solar wind B-field, OMNI)
Number of variables: 4
Number of cases/rows: 4 columns, 61 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Bx (nT)
3rd column = By (nT)
4th column = Bz (nT)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig 3B (Ion differential energy flux, Cluster-1)
Number of variables: 3
Number of cases/rows: 33 columns, 299 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Energy (keV)
Remaining 299x31 array = Ion differential energy flux (keV/cm^2/s/sr/keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig 3C (Electron differential energy flux, Cluster-1)
Number of variables: 3
Number of cases/rows: 46 columns, 892 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Energy (keV)
Remaining 892x44 array = Electron differential energy flux (keV/cm^2/s/sr/keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig 3C (Spacecraft Potential, Cluster-1)
Number of variables: 2
Number of cases/rows: 2 columns, 900 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Spacecraft potential (keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A


Fig 3D (B-field, Cluster-1)
Number of variables: 4
Number of cases/rows: 4 columns, 897 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Bx (nT)
3rd column = By (nT)
4th column = Bz (nT)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig 3E (Ion velocity, Cluster-1)
Number of variables: 4
Number of cases/rows: 4 columns, 897 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used: 
1st column = Datetime
2nd column = Vx (km/s)
3rd column = Vy (km/s)
4th column = Vz (km/s)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig 3F (Ion differential flux, Geotail)
Number of variables: 3
Number of cases/rows: 34 columns, 270 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Energy (keV)
Remaining 273x32 array = Ion differential flux (Count)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig 4A (H+ differential energy flux, Cluster-4)
Number of variables: 3
Number of cases/rows: 33 columns, 407 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used: 
1st column = Datetime
2nd column = Energy (keV)
Remaining 407x31 array = H+ differential energy flux (keV/cm^2/s/sr/keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig 4B (He+ differential energy flux, Cluster-4)
Number of variables: 3
Number of cases/rows: 33 columns, 440 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Energy (keV)
Remaining 440x31 array = He+ differential energy flux (keV/cm^2/s/sr/keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig 4C (O+ differential energy flux, Cluster-4)
Number of variables: 3
Number of cases/rows: 33 columns, 845 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Energy (keV)
Remaining 845x31 array = O+ differential energy flux (keV/cm^2/s/sr/keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig 4D (O+ by H+ differential energy flux ratio, Cluster-4)
Number of variables: 3
Number of cases/rows: 33 columns, 407 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Energy (keV)
Remaining 407x31 array = O+ by H+ differential energy flux ratio (dimensionless)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig 4E (H+ differential energy flux sorted by pitch angle, Cluster-4)
Number of variables: 3
Number of cases/rows: 18 columns, 407 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Pitch angle (deg)
Remaining 407x16 array = H+ differential energy flux (keV/cm^2/s/sr/keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig 4F (He+ differential energy flux sorted by pitch angle, Cluster-4)
Number of variables: 3
Number of cases/rows: 18 columns, 440 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Pitch angle (deg)
Remaining 440x16 array = He+ differential energy flux (keV/cm^2/s/sr/keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig 4G (O+ differential energy flux sorted by pitch angle, Cluster-4)
Number of variables: 3
Number of cases/rows: 18 columns, 845 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Pitch angle (deg)
Remaining 845x16 array = O+ differential energy flux (keV/cm^2/s/sr/keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig 4H (Electron differential energy flux sorted by pitch angle, Cluster-4)
Number of variables: 3
Number of cases/rows: 14 columns, 883 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Pitch angle (deg)
Remaining 883x12 array = Electron differential energy flux (keV/cm^2/s/sr/keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig 4I (Electron differential energy flux, 180 deg pitch angle only, Cluster-4)
Number of variables: 3
Number of cases/rows: 46 columns, 883 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Energy (keV)
Remaining 883x44 array = Electron differential energy flux (keV/cm^2/s/sr/keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig 4I (Spacecraft Potential, Cluster-4)
Number of variables: 2
Number of cases/rows: 2 columns, 900 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Spacecraft potential (keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig S1 (Walen test, Cluster 1)
Number of variables: 13
Number of cases/rows: 13 columns, 85 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = EHTx (V/m) (electric field in deHoffmann-Teller frame, X-component)
3rd column = EHTy (V/m) (electric field in deHoffmann-Teller frame, Y-component)
4th column = EHTz (V/m) (electric field in deHoffmann-Teller frame, Z-component)
5th column = ECx (V/m) (convective electric field, X-component)
6th column = ECy (V/m) (convective electric field, Y-component)
7th column = ECz (V/m) (convective electric field, Z-component)
8th column = vAx (km/s) (Alfven velocity in deHoffmann-Teller frame, X-component)
9th column = vAy (km/s) (Alfven velocity in deHoffmann-Teller frame, Y-component)
10th column = vAz (km/s) (Alfven velocity in deHoffmann-Teller frame, Z-component)
11th column = (v-vHT)x (km/s) (ion velocity in deHoffmann-Teller frame, X-component)
12th column = (v-vHT)x (km/s) (ion velocity in deHoffmann-Teller frame, Y-component)
13th column = (v-vHT)x (km/s) (ion velocity in deHoffmann-Teller frame, Z-component)
Missing data codes: N/A
Specialised formats or other abbreviations used:

Fig S2A (Electron differential energy flux, Cluster 1)
Number of variables: 3
Number of cases/rows: 46 columns, 892 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Energy (keV)
Remaining 892x44 array = Electron differential energy flux (keV/cm^2/s/sr/keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig S2A (Spacecraft Potential, Cluster-1)
Number of variables: 2
Number of cases/rows: 2 columns, 900 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Spacecraft potential (keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig S2B (Electron differential energy flux, Cluster 2)
Number of variables: 3
Number of cases/rows: 46 columns, 900 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Energy (keV)
Remaining 900x44 array = Electron differential energy flux (keV/cm^2/s/sr/keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig S2B (Spacecraft Potential, Cluster-2)
Number of variables: 2
Number of cases/rows: 2 columns, 900 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Spacecraft potential (keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig S2C (Electron differential energy flux, Cluster 3)
Number of variables: 3
Number of cases/rows: 46 columns, 899 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Energy (keV)
Remaining 899x44 array = Electron differential energy flux (keV/cm^2/s/sr/keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig S2C (Spacecraft Potential, Cluster-3)
Number of variables: 2
Number of cases/rows: 2 columns, 900 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Spacecraft potential (keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig S2D (Electron differential energy flux, Cluster 4)
Number of variables: 3
Number of cases/rows: 46 columns, 883 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Energy (keV)
Remaining 883x44 array = Electron differential energy flux (keV/cm^2/s/sr/keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig S2D (Spacecraft Potential, Cluster-4)
Number of variables: 2
Number of cases/rows: 2 columns, 900 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Spacecraft potential (keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig S3A-D (Proton differential energy flux, sunward/dawnward/duskward/antisunward, Cluster-4)
Number of variables: 3
Number of cases/rows: 33 columns, 37 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Energy (keV)
Remaining 37x31 array = Proton differential energy flux (keV/cm^2/s/sr/keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A

Fig S3E (Particle flux sum > 10 keV, sunward/dawnward/duskward/antisunward, Cluster-4)
Number of variables: 5
Number of cases/rows: 5 columns, 37 rows
Variable list, defining any abbreviations, units of measure, codes or symbols used:
1st column = Datetime
2nd column = Sunward >10 keV flux sum (keV/cm^2/s/sr/keV)
3rd column = Dawnward >10 keV flux sum (keV/cm^2/s/sr/keV)
4th column = Duskward >10 keV flux sum (keV/cm^2/s/sr/keV)
5th column = Antisunward >10 keV flux sum (keV/cm^2/s/sr/keV)
Missing data codes: N/A
Specialised formats or other abbreviations used: N/A


