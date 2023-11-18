# cop5522-project MPI stuff

Some tinkering I've done with MPI

NOTE: This code is not sorting in parallel yet.
I just wanted to gain some familiarity with MPI before I tried to modify the "samplesort" function.  
There is no comments yet - I mostly needed it clean so that I could move things around and test easier.
I've modified the "printArray" function to print to an input string instead of screen print.
I added a compareArray that will be handy later to compare a parallel sample sorted array to the "truth" serial buble sort.
Harded coded "1" in for srand() in generateRandNums - also feed it two arrays now - the second is a copy - the two will both be sorted and compared after.

I plan on making a command line argument to have some logic for srand - a number or using time(0) 
I plan on making the main function calculate the size of the array - for a CSV file input for custom problem sets - and that calculated size will be feed into the other functions rather than the command line array size.

I added "getMinMax" function - don't know if I'll use it for the final version but it gave me something to test the collectve communication MPI functions.

What the code is doing in parallel with MPI right now:
I split up the randomly generated array across the seperate processes (the last process gets the extra elements if it doesn't split evenly)
A local min and max is calculated for each process - this is reduced later to a global min and max.

I combine all the string output from all the processes and print it at the end with process 0 - this ensure everything is printed in the order that is desired - this will come in handy later.

The program time doesn't mean anything right now - I just kept it to tinker around with doing some things in process 0 versus all the processes.

Again - sorry I haven't commented the code - this is mostly to demo that I'm getting some decent proficiency with MPI.
