---
description: >-
  You will find here a compilation of the most used command while using Hydro3D.
  Your proficiency in Linux will improve greatly your CFD productivity.
---

# Tutorial

### Get information regarding a package

The following first set of commands is some of the most important to become a proficient and independent Linux user:

<table data-header-hidden><thead><tr><th></th><th></th><th data-hidden></th></tr></thead><tbody><tr><td>apt --help</td><td>Provide a concise definition of the apt package purpose and options </td><td></td></tr><tr><td>apt --man</td><td>For some package provide extensive information on the package apt</td><td></td></tr><tr><td>apt --version</td><td>Display the version of apt</td><td></td></tr></tbody></table>

### Locate the files of a package:

The command allows you to find the executable and the essential files linked to a specific package.

<table data-header-hidden><thead><tr><th></th><th></th><th data-hidden></th></tr></thead><tbody><tr><td>whereis apt</td><td>list the location of the essential files including the executable, the library, and the manual places.</td><td></td></tr><tr><td>locate apt</td><td>display all the files related to apt (extended version)</td><td></td></tr><tr><td>which apt</td><td>locate only the excutable of the package/software</td><td></td></tr></tbody></table>

### Package manager command:

It is extremely important to monitor and manage the packages installed on your workstation to avoid a forever-growing package load impeding its performance. The apt manager is intrinsic to the ubuntu distributions, but other package managers can be installed including snap, aptitude, flatpack.

<table data-header-hidden><thead><tr><th></th><th></th><th data-hidden></th></tr></thead><tbody><tr><td>sudo apt update</td><td>Update the system package and application</td><td></td></tr><tr><td>sudo apt uprade</td><td>Upgrade the version of your package and applications</td><td></td></tr><tr><td>sudo apt install</td><td>Install a given package</td><td></td></tr><tr><td>sudo apt list</td><td>List the package (dpkg) installed on your computer</td><td></td></tr><tr><td>sudo apt show</td><td>Show information regarding a given package</td><td></td></tr></tbody></table>

### Moving within Linux:

The following command allows you to move quickly through your Linux filesystem. $USER is what we call a shell variable. If your computer is named CFD1 then $USER=CFD1.&#x20;

<table data-header-hidden><thead><tr><th></th><th></th><th data-hidden></th></tr></thead><tbody><tr><td>cd </td><td>Brings you to your home directory where your main folders reside including Desktop, Document, Picture, Video, and Download.</td><td></td></tr><tr><td>cd &#x3C;path/target></td><td>Allows your to get into a specific folder.</td><td></td></tr><tr><td>cd ..</td><td>Exit a folder </td><td></td></tr><tr><td>cd /opt/</td><td>Access the /opt folder where most 3rd party apps are installed.</td><td></td></tr><tr><td>cd /media/$USER/&#x3C;target/hardrive></td><td>Access a given hard-drive</td><td></td></tr></tbody></table>

### Display content of your folder:

There is many&#x20;

```
# Get information from a program:
apt --help
apt --man
apt --version
whereis apt
locate apt
which apt

# Update & Upgrade Ubuntu:
sudo apt update
sudo upgrade
sudo apt install <Package_name>
sudo apt remove <Package_name>

# Snap update & upgrade:
sudo snap install <Package_name>
sudo snap list
sudo snap remove <Package_name>

# Moving within linux:
cd <Target_Folder>
cd ..
cd ../<Path>/<Target_Folder>
cd
cd /opt/
cd /media/$USER/<Target_Hardrive>

# Display on screen:
ls 
ll
tree 

# Add and delete files:
cp -r 
rsync -a --exclude 'tec*' <Source_Folder> <Target_Folder>
du -sh <Target_Folder>
du -sh *
rm * 
rm <target> folder 

# Combination of command:
ls <Target_Folder> | wc -l 
ps -ef | grep -rin './3dFDM.exe'
mpirun -np 10 ./3dFDM.exe | tee screen.log

# Upload and download file:
scp <target_file> user@address:path
scp user@address:path/<target_file> .

sftp user@address
sftp put -r <target_dir/file> path
sftp get -r <target_dir/file> .


# Run a simulation in ssh:
nohup mpirun -np 10 ./3dFDM.exe &
 
# Monitor your simulation:
tail -f output.dat
tail -f rms.dat
tail -f worktime.dat
htop



```

