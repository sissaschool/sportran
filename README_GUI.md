#Thermocepstrum GUI

### Installation
   See [README.md](README.md)

### Usage
   - run `thermocepstrum-gui` from the command line
   - select the input file with the  "..." button
   - select the input format. "table" is the standard column formatted data file with in the first line an header with the names of the columns. Vector data (like the 3 components of the currents) must be grouped with "c_name_[1] c_name_[2] c_name_[3]" ("c_" and "[]" must be present). Vectors will be automatically recognized.
   - "Next"
   - read the instructions and select the type of data of the various columns
   - select the units
   - "Next"
   - Input the required values. "filter width" can be left as is, it is only for visualization purposes and can be changed later
   - select the resampling frequency with the cursor. You can press "Resample" to see the effect
   - "Next"
   - the result is shown. It is possible to modify the number of cepstral coefficients and see the effect
