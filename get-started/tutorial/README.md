---
description: >-
  You will find here a compilation of the most used command while using Hydro3D.
  Your proficiency in Linux will improve greatly your CFD productivity.
---

# Tutorial

### Get information regarding a package

The following first set of commands is some of the most important to become a proficient and independent Linux user:

<table data-header-hidden><thead><tr><th></th><th></th><th data-hidden></th></tr></thead><tbody><tr><td><code>apt --help</code></td><td>Provide a concise definition of the apt package purpose and options </td><td></td></tr><tr><td><code>apt --man</code></td><td>For some package provide extensive information on the package apt</td><td></td></tr><tr><td><code>apt --version</code></td><td>Display the version of apt</td><td></td></tr></tbody></table>

### Locate the files of a package:

The command allows you to find the executable and the essential files linked to a specific package.

<table data-header-hidden><thead><tr><th></th><th></th><th data-hidden></th></tr></thead><tbody><tr><td><code>whereis apt</code></td><td>list the location of the essential files including the executable, the library, and the manual places.</td><td></td></tr><tr><td><code>locate apt</code></td><td>display all the files related to apt (extended version)</td><td></td></tr><tr><td><code>which apt</code></td><td>locate only the excutable of the package/software</td><td></td></tr></tbody></table>

### Package manager command:

It is extremely important to monitor and manage the packages installed on your workstation to avoid a forever-growing package load impeding its performance. The apt manager is intrinsic to the ubuntu distributions, but other package managers can be installed including snap, aptitude, and flatpack.

<table data-header-hidden><thead><tr><th></th><th></th><th data-hidden></th></tr></thead><tbody><tr><td><code>sudo apt update</code></td><td>Update the system package and application</td><td></td></tr><tr><td><code>sudo apt uprade</code></td><td>Upgrade the version of your package and applications</td><td></td></tr><tr><td><code>sudo apt install</code></td><td>Install a given package</td><td></td></tr><tr><td><code>sudo apt list</code></td><td>List the package (dpkg) installed on your computer</td><td></td></tr><tr><td><code>sudo apt show</code></td><td>Show information regarding a given package</td><td></td></tr></tbody></table>

### Moving within Linux:

The following command allows you to move quickly through your Linux filesystem. $USER is what we call a shell variable. If your computer is named CFD1 then $USER=CFD1.&#x20;

<table data-header-hidden><thead><tr><th></th><th></th><th data-hidden></th></tr></thead><tbody><tr><td><code>cd</code> </td><td>Brings you to your home directory where your main folders reside including Desktop, Document, Picture, Video, and Download.</td><td></td></tr><tr><td><code>cd &#x3C;path/target></code></td><td>Allows your to get into a specific folder.</td><td></td></tr><tr><td><code>cd ..</code></td><td>Exit a folder </td><td></td></tr><tr><td><code>cd /opt/</code></td><td>Access the /opt folder where most 3rd party apps are installed.</td><td></td></tr><tr><td><code>cd /media/$USER/&#x3C;TH></code></td><td>Access the targeted hard-drive</td><td></td></tr></tbody></table>

### Display the contents in your directories/files:

Displaying the content of your directories/folder correctly allows you to check the file's version, date of creation files, and permission. Here are the essentials commands:

<table data-header-hidden><thead><tr><th></th><th></th><th data-hidden></th></tr></thead><tbody><tr><td><code>ls</code></td><td>List on screen the folder content in which your are in.</td><td></td></tr><tr><td><code>ls &#x3C;TD></code></td><td>List on screen the target directory/folders</td><td></td></tr><tr><td><code>ll</code> </td><td>List on screen the detailed folder content. Size, last modified, permission.</td><td></td></tr><tr><td><code>tree &#x3C;TD></code></td><td>List in a tree format all the folder and file on screnn (sudo apt install tree).</td><td></td></tr><tr><td><code>du -sh  &#x3C;TD/TF></code></td><td>Calculate the size of folder (du -sh = disk usage -size -human).</td><td></td></tr><tr><td><code>wc -l  &#x3C;TF></code></td><td>Count and display the line count of the target file.</td><td></td></tr><tr><td><code>wc -c  &#x3C;TF></code></td><td>Count and display the character count of the target file.</td><td></td></tr><tr><td><code>cat  &#x3C;TF></code></td><td>Print the whole target file content on the screen.</td><td></td></tr><tr><td>head &#x3C;TF></td><td>Print the first 10 lines of a target file (for more option: head --help).</td><td></td></tr><tr><td><code>tail &#x3C;TF></code> </td><td>Print the last 10 lines of target file.</td><td></td></tr><tr><td><code>tail -f &#x3C;TF></code></td><td>Print the last 10 lines and print instaneously the lines if the file is updated.</td><td></td></tr></tbody></table>

### Move, Rename, Copy and delete directories/files:

The following command permits you to move, copy and delete files and directories to your convenience:

<table data-header-hidden><thead><tr><th></th><th></th><th data-hidden></th></tr></thead><tbody><tr><td><code>mv &#x3C;SF> &#x3C;TD></code></td><td>Move a source file to a target directory</td><td></td></tr><tr><td><code>mv &#x3C;file1> &#x3C;filenewname></code></td><td>Rename file1 to filenewname.</td><td></td></tr><tr><td><code>cp &#x3C;SF> ..</code></td><td>Copy the source file just outside the current directory.</td><td></td></tr><tr><td><code>cp &#x3C;file1> &#x3C;filenename></code></td><td>Copy file1 into filename</td><td></td></tr><tr><td><code>cp * &#x3C;TD></code> </td><td>Copy all the files of the current directory into the target directory. Does not copy the directory contained in the current directory.</td><td></td></tr><tr><td><code>cp -r &#x3C;SD> &#x3C;TD></code></td><td>Copy the source directory to the target directory</td><td></td></tr><tr><td><code>rsync -a --exclude 'tec*' &#x3C;SD> &#x3C;TD></code></td><td>Copy all the files of the source directory into the target directory excluding any files with the prefix 'tec' </td><td></td></tr><tr><td><code>rm *</code></td><td>Remove all the files in the current directory (Careful, the command is drastic to be used with caution)</td><td></td></tr><tr><td><code>rm &#x3C;TF></code></td><td>Remove target file</td><td></td></tr><tr><td><code>rm -r &#x3C;TD></code></td><td>Remove target directory</td><td></td></tr></tbody></table>

### Run a simulation and monitor:

<table><thead><tr><th></th><th></th><th data-hidden></th></tr></thead><tbody><tr><td><pre><code>mpirun -np 10 ./3dFDM.exe
</code></pre></td><td>Run a parallel Hydro3D simulation with 10 threads </td><td></td></tr><tr><td><pre><code>nohup mpirun -np 10 ./3dFDM.exe &#x26;
</code></pre></td><td>If you start a simulation remotely using ssh and don't want to interrupt the sim when closing the session. Start the Hydro3D simulation in the background.</td><td></td></tr><tr><td><pre><code>tail -f worktime.dat
</code></pre></td><td>Follow the time allocation to the different segments of the running simulation at each time step.</td><td></td></tr><tr><td><pre><code>tail -f rms.dat
</code></pre></td><td>Follow the mass deficit of the overall domain. </td><td></td></tr><tr><td><code>htop</code></td><td>Follow the CPU and RAM usage during the simulation (sudo apt install htop)</td><td></td></tr><tr><td><code>ps -ef</code></td><td>List the running process with the PID (Proces ID), use to kill a specific process)</td><td></td></tr><tr><td><code>kill -9 &#x3C;PID></code></td><td>Kill a specific process specific to it PID</td><td></td></tr></tbody></table>

### Work with a remote HPC:

<table data-header-hidden><thead><tr><th></th><th></th><th data-hidden></th></tr></thead><tbody><tr><td><pre><code>scp &#x3C;TF> user@address:path
</code></pre></td><td>Upload a file to a remote HPC</td><td></td></tr><tr><td><pre><code>scp -r &#x3C;TD> user@address:path
</code></pre></td><td>Upload a series of files or a directory to a remote HPC (-r=recurssive)</td><td></td></tr><tr><td><pre><code>scp user@address:path/&#x3C;TF> .
</code></pre></td><td>Download a target file from a remote HPC onto your current directory</td><td></td></tr><tr><td><pre><code>scp -r user@address:path/&#x3C;TF> .
</code></pre></td><td>Download a target folder from a remote HPC onto your current directory</td><td></td></tr></tbody></table>

### Ad:

### Ad

### Ad:

<table data-header-hidden><thead><tr><th></th><th></th><th data-hidden></th></tr></thead><tbody><tr><td><pre><code>ls &#x3C;TD> | wc -l 
</code></pre></td><td></td><td></td></tr><tr><td><pre><code>ps -ef | grep -rin './3dFDM.exe'
</code></pre></td><td></td><td></td></tr><tr><td><pre><code>mpirun -np 10 ./3dFDM.exe | tee screen.log
</code></pre></td><td></td><td></td></tr></tbody></table>

