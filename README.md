# ATLASModuleMetrology-LU
An analysis program for ATLAS ITk metrology measurements obtained from SmartScope.

Ana/ is the folder for the codes, with all of the necessary packages included. Note that HH's "date" could need installation.
  "main" is the major program, where you set the files to read, the processes to conduct and the things to output. There're many samples inside.
  "structure" sets up the SmartScope output file reading process. As the formats in Lund, Uppsala and NBI, Kopenhagen are slightly different, there's an adpation for that.
  "mathematics" has all of the necessary mathematical calculations.
    Specifically, the TwoLineIntersection calculates the middle point of the shortest line segment connecting two lines.
    MultiPointSetPlaneSVD returns a fitted plane according to the array of points.
  "measurement" is a short one with a few functions to calculate the defined position of the fifucial marks. It can be referred in my MSc thesis.
  "functions" has some of the necessary things like file line checking, reading, and the output formats.

ATLASModuleMetrology/ holds all of the outputs from the SmartScope, and the Output/ inside is for all of the outputs from this code.
