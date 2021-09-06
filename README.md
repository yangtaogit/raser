# RASER 
RAdiation SEmiconductoR 




# Build with Singularity 

Before running the code, install the Singularity on your OS. 

> ./sinularity_build.sh  

> ./singularity_run.sh 

> raser> geant4_build.sh 



# Run with Singularity 

> ./singularity_run.sh 

> raser> ./run 


# Raser unit test after you change some codes


> ./singularity_run.sh 

> raser> ./run 0.1.5
  raser> ./run 0.2.5

If the output is "Successful", the code your changed is OK. 
Otherwise, you should check the code your changed. 

