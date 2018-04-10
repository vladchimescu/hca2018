# script to run pipeline

### pancreas
loom_fn="/home/jovyan/data/tasks/how_many_cells/human_pancreas.loom"
ncell='250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,3750,4000,4250,4500,4750,5000,5500,6000,6500,7000,7500,8000,8500'
out_pfx="/home/jovyan/scratch/Group7/ashis/human_pancreas/human_pancreas_v2"
cd ~
Rscript "scratch/Group7/ashis/hca2018/pipeline.R" -loom "$loom_fn" -ncell "$ncell" -o "$out_pfx" 2>&1 | tee $out_pfx.log


### bipolar
loom_fn="/home/jovyan/data/tasks/how_many_cells/bipolar.loom"
ncell="250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,3750,4000,4250,4500,4750,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000,11000,12000,13000,14000,15000,16000,17000,18000,19000,20000,21000,22000,23000,24000,25000"
out_pfx="/home/jovyan/scratch/Group7/ashis/bipolar/bipolar_v2"
cd ~
Rscript "scratch/Group7/ashis/hca2018/pipeline.R" -loom "$loom_fn" -ncell "$ncell" -o "$out_pfx" 2>&1 | tee $out_pfx.log

