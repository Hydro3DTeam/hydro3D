# Installation

There is some prerequisite to be installed before being able to run your first simulation, and execute the tutorials:

{% tabs %}
{% tab title="Parallel XE Studio" %}
Hydro3D was developed in fortran, for academics student or staff who have access to the free version of Intel Parallel XE studio it is advice to install the fortran and mpi compiler from the free version:

In the following, we install the Intel Parallel Studio XE 2018 4th version as it is the version currently used by the whole research group but feel free to try out newer version and let us know if it is successful in compiling Hydro3D.

Click [Parallel Studio XE](https://software.intel.com/content/www/us/en/develop/tools/parallel-studio-xe/choose-download/student-linux.html) to register and download your version and its academic code. Once it is downloaded in your Download folder. Proceed to the following command:

```
sudo apt-get install g++
```

Decompressed the Parallel Studio XE folder:

```bash
tar -xzf parallel_studio_xe_2018_update4_cluster_edition.tgz
```

Go in the decompressed the folder and started the GUI installation:

```bash
cd parallel_studio_xe_2018_update4_cluster_edition
sudo ./install_GUI.sh
```

![](<../.gitbook/assets/Screenshot from 2020-11-10 16-37-11.png>)

![](<../.gitbook/assets/Screenshot from 2020-11-10 16-35-49.png>)

![](<../.gitbook/assets/Screenshot from 2020-11-10 16-38-36.png>)

![](<../.gitbook/assets/Screenshot from 2020-11-10 16-39-09.png>)

![](<../.gitbook/assets/Screenshot from 2020-11-10 16-40-16.png>)

![](<../.gitbook/assets/Screenshot from 2020-11-19 09-49-18.png>)

Congrats, you have just installed the compiler ! Now if we want to use in terminal to compile Hydro3D, it is important to do some change in your .bashrc. Do not worry for beginners in Linux it might looks a bit scary at first. The .bashrc is the file that allows you to customize your terminal experience such as adding alias to shorcut some of the actions you want execute. To access the .bashrc enter the following commands with your preferred text editor here I use nano:

```
cd ~
nano .bashrc
```

In order to access the Intel Fortran and MPI libraries from the terminal, it is necessary to add these lines at the end of your .bashrc file.

```
source /opt/intel/bin/compilervars.sh intel64                              # setup the mpi compiler environment
source /opt/intel/bin/ifortvars.sh intel64                                 # setup the fortran compiler environment
export PATH=$PATH:/opt/intel/compilers_and_libraries/linux/mpi/intel64/bin # setup the path 
export I_MPI_F90=ifort                                                     # specify the fortran compiler to use for compiling parallel code 
```

Once the .bashrc is saved and close in the terminal enter to update apply the change of .bashrc within across the terminal:

```
cd ~
. .bashrc
reset
```

\*\* Optional Section Start\
Prior to testing that the compilers and environment are correctly set up, for the linux enthusiam I will give you a bonus optional task. Reopen your .bashrc and copy the following:

```
# Open .bashrc with nano from any location you might find your 
alias bashrc='nano ~/.bashrc'  

# Update the bashrc and reset the terminal window (';' allows to enter two commands):
alias bashsync='. ~/.bashrc; reset'     
```

Update for the last time your .bashrc with the normal command. What we have done is create alias to enter your .bashrc regarding where your are located in your computer ( just type bashrc) and to update it (bashsync).\
Optional Section End \*\*

Before running the Hydro3D tutorials, let's make sure the libraries can be called from your terminal, for that please try to run the two little programs. For the fortran library run this program:

```
ifort parallecode.for -o parallecode.exe
mpirun -np 4 parallelcode.exe
```


{% endtab %}
{% endtabs %}
