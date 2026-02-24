ListingAnalyst Version 1.3

This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the software.
 
DOCUMENTATION:
Winston, R. B., and Paulinski, Scott. 2013. ListingAnalyst: A Program for Analyzing the Main Output File from MODFLOW. Ground Water. doi: 10.1111/gwat.12054.

INSTALLATION:
The program is distributed as a .zip file. The program can be installed by extracting it from the .zip file.

RUNNING:
To run ListingAnalyst, simply double-click on ListingAnalyst.exe from the bin directory of the distribution file.

COMPILING
ListingAnalyst is compiled with Delphi 11.3 from Embarcadero. 
http://www.embarcadero.com/ 

ListingAnalyst uses a number of custom components that must be installed 
in Delphi 11.3 before compiling ListingAnalyst.  Some are included with the 
ModelMuse source code.  Additional required files or components are 
listed below.  

General instructions for installing packages in Delphi 11.3
1. If the component comes with an installer run the installer.
2. If you are compiling the components from source code, you need to add the 
directories containing the source code to the Library path for both the 
Windows 32 and Windows 64 bit platforms. (Tools|Options) 
then look in "Environment Options|Delphi Options|Library".
3. If you are compiling the components from source code and the components 
are separated into run-time and design-time packages, you build the runtime 
package and then build and install the design-time package.

Install JCL and JVCL They can be obtained from http://www.delphi-jedi.org/
Add the following JCL directories to the Library path if they are not added
automatically when installing the JCL.
source\common
source\windows

Get and install VirtualTreeView version 5.0.0 or later.
http://www.soft-gems.net/

