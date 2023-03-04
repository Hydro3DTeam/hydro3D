# Tutorial

```
# Get information from program:
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

# Display:
ls 
ll





cp -r 
rsync -a --exclude 'tec*' <Source_Folder> <Target_Folder>
du -sh <Target_Folder>
du -sh *
rm * 
rm <target> folder 


ls <Target_Folder> | wc -l 
ps -ef | grep -rin './3dFDM.exe'
mpirun -np 10 ./3dFDM.exe | tee screen.log

```

