### This script is meant to perform virtual screening of compounds against TREX structures, using a standard .smi input file for the compounds.
### Set v= to the name of your input .smi file, which should be in the 'AutoDock' directory.
v='TIM218'
### Set g= to the number of cores desired for the screen.
let g=3
### Set modes= to the number of binding poses you want reported
let modes=9
### Select one of the 3 following TREX structures to dock against, partially or completely, by defining t=.
# 1 = 'trex1-dAMP', active-site
# 2 = 'trex2-apo', active-site
# 3 = 'trex1-apo', active-site
# 4 = 'trex1-dAMP', whole structure
# 5 = 'trex2-apo', whole structure
# 6 = 'trex1-apo', whole structure
# 7 = 'mT1-ssDNA', whole structure
# 8 = 'mT1-Mg', whole structure
# 9 = 'hT1-Mg', whole structure
# 10 = 'RAD54', whole structure
# 11 = 'RAD54', DNA-binding site
# 12 = 'hT1-Mg', active-site
# 13 = 'hT1', active-site
# 14 = 'hT1', whole structure
let t=12
### Set p= to the desired exhaustiveness for the screen (standard =8).
let p=100

#####

cd ~/Documents/Research/Perrino\ Lab/AutoDock/
mkdir "$v"_Screen
mkdir "$v"_Screen/vinalogs
mkdir "$v"_Screen/vinadocks
mkdir "$v"_Screen/HITS
cd "$v"_Screen/
cp ../Scripts/screen.sh script.sh
mv ../"$v".smi .
x=$(cat "$v".smi | wc -l)
for((i=1;i<="$x";i++));
do
s=$(awk -v u="$i" '{if (NR==u){print $2}}' "$v".smi)
babel -i smi "$v".smi -o pdb "$v".pdb  -f $i -l $i --gen3d
/Library/MGLTools/1.5.6/bin/pythonsh /Library/MGLTools/1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l "$v".pdb -o "$v".pdbqt
if [ "$t" -eq 1 ]
then
    vina --receptor ~/Documents/Research/Perrino\ Lab/AutoDock/trex1-dAMP.pdbqt --ligand "$v".pdbqt --center_x -10 --center_y -5 --center_z 5 --size_x 25 --size_y 25 --size_z 25 --cpu "$g" --num_modes "$modes" --exhaustiveness "$p" --log vinalogs/log_"$v"_"$i"_"$s" --out vinadocks/"$v"_"$i"_"$s"_out.pdbqt
fi
if [ "$t" -eq 2 ]
then
    vina --receptor ~/Documents/Research/Perrino\ Lab/AutoDock/trex2-apo.pdbqt --ligand "$v".pdbqt --center_x -10 --center_y -5 --center_z 5 --size_x 25 --size_y 25 --size_z 25 --cpu "$g" --num_modes "$modes" --exhaustiveness "$p" --log vinalogs/log_"$v"_"$i"_"$s" --out vinadocks/"$v"_"$i"_"$s"_out.pdbqt
fi
if [ "$t" -eq 3 ]
then
    vina --receptor ~/Documents/Research/Perrino\ Lab/AutoDock/trex1-apo.pdbqt --ligand "$v".pdbqt --center_x -10 --center_y -5 --center_z 5 --size_x 25 --size_y 25 --size_z 25 --cpu "$g" --num_modes "$modes" --exhaustiveness "$p" --log vinalogs/log_"$v"_"$i"_"$s" --out vinadocks/"$v"_"$i"_"$s"_out.pdbqt
fi
if [ "$t" -eq 4 ]
then
    vina --receptor ~/Documents/Research/Perrino\ Lab/AutoDock/trex1-dAMP.pdbqt --ligand "$v".pdbqt --center_x 3 --center_y 2 --center_z 21 --size_x 75 --size_y 75 --size_z 75 --cpu "$g" --num_modes "$modes" --exhaustiveness "$p" --log vinalogs/log_"$v"_"$i"_"$s" --out vinadocks/"$v"_"$i"_"$s"_out.pdbqt
fi
if [ "$t" -eq 5 ]
then
    vina --receptor ~/Documents/Research/Perrino\ Lab/AutoDock/trex2-apo.pdbqt --ligand "$v".pdbqt --center_x 3 --center_y 2 --center_z 21 --size_x 75 --size_y 75 --size_z 75 --cpu "$g" --num_modes "$modes" --exhaustiveness "$p" --log vinalogs/log_"$v"_"$i"_"$s" --out vinadocks/"$v"_"$i"_"$s"_out.pdbqt
fi
if [ "$t" -eq 6 ]
then
    vina --receptor ~/Documents/Research/Perrino\ Lab/AutoDock/trex1-apo.pdbqt --ligand "$v".pdbqt --center_x 3 --center_y 2 --center_z 21 --size_x 75 --size_y 75 --size_z 75 --cpu "$g" --num_modes "$modes" --exhaustiveness "$p" --log vinalogs/log_"$v"_"$i"_"$s" --out vinadocks/"$v"_"$i"_"$s"_out.pdbqt
fi
if [ "$t" -eq 7 ]
then
    vina --receptor ~/Documents/Research/Perrino\ Lab/AutoDock/mT1-ssDNA.pdbqt --ligand "$v".pdbqt --center_x 20 --center_y 10 --center_z -20 --size_x 75 --size_y 75 --size_z 75 --cpu "$g" --num_modes "$modes" --exhaustiveness "$p" --log vinalogs/log_"$v"_"$i"_"$s" --out vinadocks/"$v"_"$i"_"$s"_out.pdbqt
fi
if [ "$t" -eq 8 ]
then
vina --receptor ~/Documents/Research/Perrino\ Lab/AutoDock/mT1-Mg.pdbqt --ligand "$v".pdbqt --center_x 3 --center_y 2 --center_z 21 --size_x 75 --size_y 75 --size_z 75 --cpu "$g" --num_modes "$modes" --exhaustiveness "$p" --log vinalogs/log_"$v"_"$i"_"$s" --out vinadocks/"$v"_"$i"_"$s"_out.pdbqt
fi
if [ "$t" -eq 9 ]
then
vina --receptor ~/Documents/Research/Perrino\ Lab/AutoDock/hT1-Mg.pdbqt --ligand "$v".pdbqt --center_x 3 --center_y 2 --center_z 21 --size_x 75 --size_y 75 --size_z 75 --cpu "$g" --num_modes "$modes" --exhaustiveness "$p" --log vinalogs/log_"$v"_"$i"_"$s" --out vinadocks/"$v"_"$i"_"$s"_out.pdbqt
fi
if [ "$t" -eq 10 ]
then
vina --receptor ~/Documents/Research/Perrino\ Lab/AutoDock/RAD54_ed1.pdbqt --ligand "$v".pdbqt --center_x 23 --center_y 29 --center_z 22 --size_x 82 --size_y 100 --size_z 72 --cpu "$g" --num_modes "$modes" --exhaustiveness "$p" --log vinalogs/log_"$v"_"$i"_"$s" --out vinadocks/"$v"_"$i"_"$s"_out.pdbqt
fi
if [ "$t" -eq 11 ]
then
vina --receptor ~/Documents/Research/Perrino\ Lab/AutoDock/RAD54_ed1.pdbqt --ligand "$v".pdbqt --center_x 26 --center_y 50 --center_z 32 --size_x 38 --size_y 39 --size_z 35 --cpu "$g" --num_modes "$modes" --exhaustiveness "$p" --log vinalogs/log_"$v"_"$i"_"$s" --out vinadocks/"$v"_"$i"_"$s"_out.pdbqt
fi
if [ "$t" -eq 12 ]
then
vina --receptor ~/Documents/Research/Perrino\ Lab/AutoDock/hT1-Mg.pdbqt --ligand "$v".pdbqt --center_x -10 --center_y -5 --center_z 5 --size_x 25 --size_y 25 --size_z 25 --cpu "$g" --num_modes "$modes" --exhaustiveness "$p" --log vinalogs/log_"$v"_"$i"_"$s" --out vinadocks/"$v"_"$i"_"$s"_out.pdbqt
fi
if [ "$t" -eq 13 ]
then
vina --receptor ~/Documents/Research/Perrino\ Lab/AutoDock/hT1.pdbqt --ligand "$v".pdbqt --center_x -10 --center_y -5 --center_z 5 --size_x 25 --size_y 25 --size_z 25 --cpu "$g" --num_modes "$modes" --exhaustiveness "$p" --log vinalogs/log_"$v"_"$i"_"$s" --out vinadocks/"$v"_"$i"_"$s"_out.pdbqt
fi
if [ "$t" -eq 14 ]
then
vina --receptor ~/Documents/Research/Perrino\ Lab/AutoDock/hT1.pdbqt --ligand "$v".pdbqt --center_x 3 --center_y 2 --center_z 21 --size_x 75 --size_y 75 --size_z 75 --cpu "$g" --num_modes "$modes" --exhaustiveness "$p" --log vinalogs/log_"$v"_"$i"_"$s" --out vinadocks/"$v"_"$i"_"$s"_out.pdbqt
fi
awk '{if (NF==4 && $1==1){print $2,FILENAME}}' vinalogs/log_"$v"_"$i"_"$s" >> scores.dat
awk -v q="vinadocks/"$v"_"$i"_"$s"_out.pdbqt" '{if (NF==4 && $2<-9.9 && $1==1){print $2,q}}' vinalogs/log_"$v"_"$i"_"$s" >> hits.dat
done
n=$(cat hits.dat | wc -l)
for ((i=1;i<="$n";i++));
do
temp=$(awk -v m="$i" '{if (NR==m){print $2}}' hits.dat)
cp "$temp" HITS/
done
mv hits.dat HITS/hits.dat
rm "$v".pdb
rm "$v".pdbqt
cd ~
