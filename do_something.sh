#/bin/bash


# for ijob in $(qstat -u petavazohi | grep MBJ | awk '{ print $1 }' | awk -F '.' '{ print $1}')
# do
#     qdel $ijob
# done

for idir in mp-*
do
    rm ${idir}/MBJ/*
done




# mkdir MoS2-bulk
# for idir in [0-9].*
# do
#     echo $idir
#     cp ${idir}/bands/POSCAR MoS2-bulk/${idir}.POSCAR
# done



# mkdir /tmp/MoS2
# for idir in 0 1 11 12 13 14  15  2  3  4  5  6  7  8  9
# do
#     echo $idir
#     mkdir /tmp/MoS2/${idir}
#     cp $idir/bands/CONTCAR /tmp/MoS2/${idir}
# done

    









# for idir in LDA  PBE  PBEsol
# do
#     for jdir in mp-13_Fe  mp-19035_BaFeO3  mp-21078_Fe3Ge  mp-510624_SrFeO3  mp-778_Fe2P
#     do
# 	for kdir in DFT+U  regular
# 	do
# 	    for ldir in bands  dos  relax  scf
# 	    do
# 		cp MCMC_selected/${idir}/${jdir}-${idir}/${kdir}/${ldir}/CONTCAR DataAvailnpj/exploration_stage/${idir}/${jdir}/${kdir}/${ldir}/CONTCAR
# 		echo ${idir}-${jdir}-${kdir}-${ldir}
# 	    done
# 	done
#     done
# done

# for idir in LDA  PBE  PBEsol
# do
#     for jdir in mp-1105372-Fe5SiB2   mp-3805-AlFe2B2   mp-9913-Fe5PB2
#     do
# 	for kdir in DFT+U  regular
# 	do
# 	    cp evaluation_stage/${jdir}/${idir}/${kdir}/CONTCAR DataAvailnpj/1evaluation_stage/${idir}/${jdir}/${kdir}/relax/CONTCAR
# 	    echo ${idir}-${jdir}-${kdir}-${ldir}
# 	done
#     done
# done


# for idir in LDA  PBE  PBEsol
# do
#     for jdir in mp-568961-BaFe2As2
#     do
# 	for kdir in DFT+U   regular 
# 	do
# 	    for ldir in AFM-checkerboard  AFM-striped  FM
# 	    do
# 		cp  evaluation_stage/${jdir}/${idir}/${kdir}/${ldir}/INCAR  DataAvailnpj/evaluation_stage/${idir}/${jdir}/${kdir}/${ldir}/INCAR
# 		echo ${idir}-${jdir}-${kdir}-${ldir}
# 	    done
# 	done
#     done
# done


# For idir in LDA PBE PBEsol
# do
#     for jdir in DFT+U  regular
#     do
# 	for kdir in AFM-checkerboard  AFM-striped  FM
# 	do
# 	    echo ${idir}-${jdir}-${kdir}
# 	    grep TOTEN ${idir}/${jdir}/${kdir}/OUTCAR | tail -1 | awk '{print $5}'
# 	done
#     done
# done




# for idir in LDA #PBE PBEsol
# do
#     cd ${idir}
#     for jdir in DFT+U  #regular
#     do
# 	cd ${jdir}
# 	for kdir in AFM-checkerboard  AFM-striped  FM
# 	do
# 	    cd ${kdir}
# 	    echo ${idir}-${jdir}-${kdir}
# 	    get_calculation_info.py > ${idir}-${jdir}-${kdir}.log
# 	    mv ${idir}-${jdir}-${kdir}.log ../../../
# 	    cd ..
# 	done
# 	cd ..
#     done
#     cd ..
# done

	    
# for ipot in $(find . -name POTCAR)
# do
#     grep 'TITEL\s*=\s*PAW' $ipot >> ${ipot}_options
# done

   

# for idir in LDA PBE PBEsol
# do
#     for jdir in mp-1105372-Fe5SiB2 mp-18905-FeO mp-19770-Fe2O3 mp-3805-AlFe2B2 mp-568961-BaFe2As2 mp-9913-Fe5PB2
#     do
# 	for kdir in LDA PBE PBEsol
# 	do
# 	    if [ $jdir == 'mp-19770-Fe2O3' ]
# 	    then
# 		if [ ${kdir} == "LDA" ]
# 		then
# 		    rm ${idir}/${jdir}/PBE*
# 		    mv ${idir}/${jdir}/LDA+U1/  ${idir}/${jdir}/DFT+U/
# 		    mv ${idir}/${jdir}/LDANoU/  ${idir}/${jdir}/regular
# 		elif [ ${kdir} == "PBE" ]
# 		then
# 		    rm ${idir}/${jdir}/LDA*
# 		    rm ${idir}/${jdir}/PBEsol*
# 		    mv ${idir}/${jdir}/PBE+U1/  ${idir}/${jdir}/DFT+U/
# 		    mv ${idir}/${jdir}/PBENoU1/  ${idir}/${jdir}/regular
# 		else
# 		    rm ${idir}/${jdir}/LDA*
# 		    rm ${idir}/${jdir}/PBE+*
# 		    rm ${idir}/${jdir}/PBENo*
# 		    mv ${idir}/${jdir}/PBEsol+U1/  ${idir}/${jdir}/DFT+U/
# 		    mv ${idir}/${jdir}/PBEsolNoU1/  ${idir}/${jdir}/regular
# 		fi
# 	    else
# 		if [ $idir != $kdir ]
# 		then
# 		    # rm -r ${idir}/${jdir}/${kdir}
# 		    # mv ${idir}/${jdir}/${idir}/* ${idir}/${jdir}/
# 		    # rm -r ${idir}/${jdir}/${idir}
# 		    echo " "
# 		fi
# 	    fi
# 	done
#     done
# done

