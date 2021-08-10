#/bin/bash

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

