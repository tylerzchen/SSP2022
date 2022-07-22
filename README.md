------------------ Orbital Determination Code, SSP 2022 ------------------
Name: Tyler Chen
Date: July 21st, 2022
Asteroid: 1999 GJ2
Team: Team 01 - Glacier National Park
Teammates: Anoushka Chitnis, Sam Kleiman, Alison Soong
--------------------------------------------------------------------------
API
--------------------------------------------------------------------------
Generated using https://pypi.org/project/pdoc/

Respective files for the API Documentation are with the ChenOD folder, however, the HTML file can't be opened independently if the local host is not running. Due to this you must follow the following instructions.

Instructions:
    - Install pdoc using pip install in the terminal:
    
        pip install pdoc
    
    - Download ChenOD package
    - Open terminal
    - Navigate to ChenOD directory
    - Type the following in the terminal:
    
          pdoc --http localhost:XXXX ODLIB.py
          
          where XXXX is any 4 digit number
          
    - If your 4 digit number does not work, try another 4 digit number

Once your local host is running you can update the API yourself with your own edits to the ODLIB.py docstrings and generate your own HTML files. 
You can also open the HTML file included, but this can only be done with the 2344 local host. 
--------------------------------------------------------------------------
OD Generation
--------------------------------------------------------------------------
To perform orbital determination, ephemeris generation, and/or monte carlo simulation with your own inputs, you can use the ODGeneration.ipynb file from the ChenOD folder.

File Formatting:

To ensure proper orbital determination generation your data input and sun input files must be formatted correctly. All observation times should be in julian days. All right acensions should be in HMS. All declinations should be in DMS. Sun vectors should be in AU, equatorial plane, J200, apparent from JPL.

Your data input .txt file from your three observations for orbital determination should be formatted as follows:

Obs1 Time (Julian Days), RA1 Hours, RA1 Minutes, RA1 Seconds, Dec1 Degrees, Dec1 Arcminutes, Dec1 Arcseconds, R1X, R1Y, R1Z
Obs2 Time (Julian Days), RA2 Hours, RA2 Minutes, RA2 Seconds, Dec2 Degrees, Dec2 Arcminutes, Dec2 Arcseconds, R2X, R2Y, R2Z
Obs3 Time (Julian Days), RA3 Hours, RA3 Minutes, RA3 Seconds, Dec3 Degrees, Dec3 Arcminutes, Dec3 Arcseconds, R3X, R3Y, R3Z

Ex.

2459758.6900936, 16, 29, 13.23, 11, 49, 50.5, -1.094708696204573E-01,  9.273413097249199E-01,  4.019595794929602E-01
2459772.6782503, 16, 22, 53.90, 11, 23, 23.8, -3.396653695756806E-01,  8.791228948267975E-01,  3.810538933651598E-01
2459774.6955906, 16, 22, 48.54, 11, 10, 0.9, -3.716017622155263E-01,  8.681063517238574E-01,  3.762768881806196E-01

Your suns input .txt file for ephemeris generation can contain as many times and sun vectors as desired and should be formatted as follows:

Obs1 Time (Julian Days), R1X, R1Y, R1Z
Obs2 Time (Julian Days), R2X, R2Y, R2Z
.
.
.
.
ObsN Time (Julian Days), RNX, RNY, RNZ

Ex.

2459758.6900936, -1.094708696204573E-01,  9.273413097249199E-01,  4.019595794929602E-01
2459772.6782503, -3.396653695756806E-01,  8.791228948267975E-01,  3.810538933651598E-01
2459774.6955906, -3.716017622155263E-01,  8.681063517238574E-01,  3.762768881806196E-01

After creating your dataInput.txt file and sunsInput.txt file in the proper formatting, follow the instructions in the ODGeneration.ipynb file for OD generation.
--------------------------------------------------------------------------
Asteroid Orbit Visualization
--------------------------------------------------------------------------   
Visualization generated with https://docs.poliastro.space/en/stable/ and https://www.astropy.org/

In order to run visualization software you must have the correct versions of astropy and poliastro.
If you don't have the astropy or poliastro packages installed, follow the following instructions.

Instructions:
    - Install the *CORRECT* version of astropy using pip install in the terminal:
    
        pip install astropy==4.0
    
    - Install the *CORRECT* version of poliastro using pip install in the terminal:

        pip install poliastro==0.14.0

After installing the correct versions of astropy and poliastro, follow the directions inside the second kernel of the ODGeneration.ipynb file from the ChenOD folder.

If you already have astropy or poliastro installed, check to make sure you have the correct version using:

    pip list

If you have the incorrect versions of astropy and poliastro you will have to do the following in the terminal:

    pip uninstall astropy
    
    pip uninstall poliastro

After you have uninstalled your incorrect versions of astropy and poliastro, you can reinstall both packages using the instructions for installing the correct versions.
--------------------------------------------------------------------------
