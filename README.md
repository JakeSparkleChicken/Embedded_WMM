# Embedded_WMM
Port of the World Magnetic Model code to embedded microcontrollers

The WMM program is used to determine the exact magnetic declination for any spot on the planet.  Unfortunately, it is a really beefy bit of code if you're trying to mash it into a microcontroller.  The geoid information alone is 4MB of data.  Try fitting that into a Cortex-M!
The first thing dispensed with was the geoid code.  I then excised all of the functions related to asking users for input and displaying the information on a terminal.  The last thing that I did was move the coefficients from a file to a function.  Unfortunately, that pushes the size of the executable above what an STM32F030R8 can handle.  I'll be doing the rest of the development on an STM32F4 or STM32L4.
