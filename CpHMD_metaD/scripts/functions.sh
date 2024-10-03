# This file contains all the functions called by CpH-MD.sh script
#It also contains the definition of some constants.
#
# Please don't delete it
#
#######################################################################
# SC 28-11-2011
# In this part some constants can be defined which may change with 
# gromacs version. When adapting to new gromacs versions, please check 
# in here.
#
# Needed to know how to delimit moleculetype section in pdb2gmx output 
# topology. 
# Careful: it may change with gromacs version
# This one works with 3.2.1 and 4.0.7 at least
export TOPSTRING="; Include water topology" 
#
#######################################################################
#
correct_variables ()
{                   
    offset=1000 # which bypasses the user defined value (older settings)
    GridSize=`awk -v i=$GridSize 'BEGIN{print (i==0) ? 0 : i}'`
    offset=`awk -v i=$offset 'BEGIN{print (i==0) ? 0 : i}'` 
}

build_forcefield ()
{
    # Call the modified forcefield:
    ln -s -f "$CpHDIR"/top/${ffID}.ff .; ln -s -f "$CpHDIR"/top/*.dat .
}

check_files ()
{
    # Check for other necessary files (filenames must be respected)
    for f in ${runname}.mdp $GROin $TOPin $GroDIR $PDBin $NDXin
      do
      if [ ! -f $f ]; then
          echo message E  "File $f is missing!!!... Program will crash"
      fi
    done
    #Check topology file for delimiting strings ;SC 28-11-2011
    awk -v p=$prog '
    BEGIN{
    b_string=";[ begin_tit_molecule ]"; 
    e_string=";[ end_tit_molecule ]";
    b=0;
    end=0;
    order=1;
    warnfmt=(p": Warning: String %s appears repeated in %s file. First occurrence will be considered.\n");
    errfmt=(p": Error: String %s is missing in %s file.\n")
    }
    /^ *; *\[ +begin_tit_molecule +\]/{
    if (b) printf warnfmt, b_string, ARGV[1]; else b=1}
    /^ *; *\[ +end_tit_molecule +\]/{
    if (end) printf warnfmt, e_string, ARGV[1]; 
    else {if (!b) order=0; end=1}}
    END{
    if (!b) {printf errfmt, b_string, ARGV[1]; exit 1};
    if (!end) {printf errfmt, e_string, ARGV[1]; exit 1};
    if (!order) {printf (p": Error: String %s is before string %s in %s file.\n"), e_string, b_string, ARGV[1]; exit 1}
    }' $TOPin >&2  
    #
    # Check directories
    for d in $PetitDIR $StDIR $DelphiDir
      do
      if [ ! -d $d ]; then 
          message E  "Directory \$d does not exist!!!... Program will crash"
      fi
    done    
    if [[ ! -f $RULEin ]]; then
        message E  "File \$RULEin is missing!!!... Program will crash"
    fi

    if [[ $PBdim -ne 0 && $PBdim -ne 2 ]]; then
        message E  "Number of dimensions ($PBdim) is not supported in ${CpHModule} module!!!... Program will crash"
    fi
    #
    # Detect cutoff-scheme to be used
    scheme=`awk '/cutoff-scheme/ {print $3}' ${runname}.mdp `
    if [[ $scheme == "verlet" ]]; then
	message W "Scheme used is verlet."
    elif [[ $scheme == "group" ]]; then
        message W "Scheme used is Group based."
	## Detect if the gromacs version supports group based
	#
	$GroDIR help  > ./gro-ver.txt 2>&1 
	version=`awk '/GROMACS:/ {print $(NF)}' gro-ver.txt`
	if [[ $version != "5.1.5" ]]; then
	    message E  "Cut-off scheme ${scheme} is not supported in ${GroDIR} module!!!... Program will crash"
	fi	
    else
        message W "Scheme is not group nor verlet... Performing method as group based cut off"
    fi
    ###########################################################################
    ## Check if gromacs version is not 2020.2 since it has freezedims bugged ##
    ###########################################################################
    $GroDIR help  > ./gro-ver.txt 2>&1 
    version=`awk '/GROMACS:/ {print $(NF)}' gro-ver.txt`
    if [[ $version == "2020.2" ]]; then
	message E  "Gromacs ${version} is not supported for CpHMD since freezedims is bugged !!... Program will crash"
    fi	
    ### Adition 11/05/2023 Nolive ###
    ## Check if all includes in the topology have been copied to the folder
    #
    for include in `awk '/^ *; *\[ +begin_tit_molecule +\]/{exit}; /#include/ {print $2}'  $TOPin ` `awk 'n;/^ *; *\[ +end_tit_molecule +\]/{n=1}' $TOPin | awk '/#include/ {print $2}' ` ;do
	if [ !  -f `echo $include | sed 's/"//g'` ]; then
	    message E  "$include file not found... Program will crash"
	fi
    done
    #####################################
    ## Check if plumed is being used ##
    #####################################
    # To Be updated to check for multiple versions
    if [[ $plumed != "no" ]] && [[ ! -f ${runname}_plumed.dat ]] ;then
    message E  "File ${runname}_plumed.dat not found... Program will crash"        
    fi


}

make_auxiliary_files ()
{
    # Rename original files
    cp -f $TOPin TMP_CpHMD.top
    cp -f $GROin TMP_effective.gro
    cp -f $NDXin TMP_CpHMD.ndx
    #
    # Create auxiliary files with the unchanged part of the topology file; SC 28-11-2011
    awk '/^ *; *\[ +begin_tit_molecule +\]/{exit}{print}' $TOPin > TMP_CpHMD.topbegin 
    awk 'n;/^ *; *\[ +end_tit_molecule +\]/{n=1}' $TOPin > TMP_CpHMD.topend
  
    # Check if there are different ionic strength values in MD and PB
    if [[ $ionicstrMD == "" ]]; then ionicstrMD=$ionicstr; fi
    # Convert Ionic Strength from Molar to molecule/nm^3.
    ionicstrMolecule=$(awk -v i=${ionicstrMD} 'BEGIN{print i*0.6022}')
    #
    # Correct Ionic Strength and Temperature in the .mdp file according 
    # to your parameters
    sed "s/\(ionicstrength *= \).*\$/\1 $ionicstrMolecule/" \
        ${runname}.mdp > TMP_aux1.mdp
    message W "Ionic strength in ${runname}.mdp file was changed to $ionicstrMolecule molecules/nm^3 (which corresponds to $ionicstr Molar)."
    #
    tcgroups=`awk '$1=="tc-grps"{for (i=3;i<=NF;i++) if ($i !~ /^;/) {n++; if ($i ~ /;/) exit} else exit}END{print n}' ${runname}.mdp` #SC 25-11-2011
    awk -v e=$tcgroups -v t=$temp '!/^ref_t *=/{print $0}; 
        /^ref_t *=/{printf "ref_t               =  "; 
        for (i=1;i<=e;i++) printf "%-9s", t ;print ""}' \
            TMP_aux1.mdp > ${runname}.mdp #SC 25-11-2011
        message W "Temperature in ${runname}.mdp file was changed to $temp."

    # Make relaxation .mdp gromacs file
    # Remove constraints & change NPT to NVT
    sed "s/\(nsteps *= \).*\$/\1 $RelaxSteps/ 
         s/\(Pcoupl *= \).*\$/\1 No/" \
        ${runname}.mdp > TMP_relax.mdp
    #
    #### Adition on 13/01/2021 to run in either group and verlet cutoff based runs ####
    ### Add Different protocol depending on the cut-off scheme ###
    # Check which cut-off is being used #
    scheme=`awk '/cutoff-scheme/ {print $3}' ${runname}.mdp `
    if [[ $scheme == "verlet" ]]
    then
        echo -e "\nfreezegrps          =  Solute\nfreezedim           =\
        Y Y Y\n" >> TMP_relax.mdp
    elif [[ $scheme == "group" ]]
    then
        echo -e "\nfreezegrps          =  Solute\nfreezedim           =\
        Y Y Y\n\nenergygrp_excl      =  Solute Solute" >> TMP_relax.mdp
    else
        message W "Scheme not group nor verlet... Performing method as group based cut off"
        echo -e "\nfreezegrps          =  Solute\nfreezedim           =\
        Y Y Y\n\nenergygrp_excl      =  Solute Solute" >> TMP_relax.mdp

    fi
    #
    # Make effective .mdp gromacs file
    sed  "s/\(nsteps *= \).*\$/\1 $EffectiveSteps/" ${runname}.mdp > TMP_effective.mdp
    #
    # Modify Effective MDP file in case of using Position Restraints
    if [[ -f $PosRe && $PosRe != *(" ") ]]; then 
        sed -i "s/\(define *= \).*\$/\1 -DPOSRES/" TMP_effective.mdp
    else
        if [ ! -f $PosRe ]; then message W \
        "File $PosRe is missing. Position Restraints will not be used"; fi
    fi
    #
    ### 24/04 changed -p $TOPin to -p TMP_CpHMD.top for the grompp
    ##to run using the folder's FF
    "$GroDIR" grompp -f TMP_effective.mdp -po TMP_effective_out.mdp \
        -c  $PDBin -p TMP_CpHMD.top -pp TMP_processed.top -n TMP_CpHMD.ndx \
        -o TMP_CpHMD.tpr -maxwarn 1000 -quiet
    if [ $? -eq 1 ]; then
        echo "Job cancelled due to process error."
        scancel ${SLURM_JOB_ID}
    fi
    rm -f TMP_aux*.mdp
    #
 
    echo "Protein" | "$GroDIR" editconf -f $PDBin \
			       -o TMP_protein.pdb -n TMP_CpHMD.ndx -quiet
    
    echo "Solute" | "$GroDIR" editconf -f $PDBin \
			       -o TMP_solute.pdb -n TMP_CpHMD.ndx -quiet
    ##################################################################
    ### Correct the Cter O1 atom name so it doesn't crash later on ###
    ##################################################################
    Ctres=`awk 'BEGIN{chain_old=" "};/ATOM/{if ($3 == "C") last_c=substr($0,23,4); chain_now=substr($0,22,1); if (chain_now != chain_old) {if (chain_old != " ") ; chain_old=chain_now}}END{if (chain_old != " ") print last_c+0}' TMP_protein.pdb`
    sed -i "/ $Ctres /s/ O1 / O /" TMP_protein.pdb
    #################################################################
    ### Getting the correct pbp file for PBMC                     ###
    #################################################################
    if [[ $PBdim == 0 ]]
    then
	sed "s/\(export ionicstr*=\).*\$/\1$ionicstr/
             s/\(export epsin*=\).*\$/\1$epsin/" \
		 "$DelphiDir"/DELPHI_Prot.pbp > DELPHI.pbp
    elif [[ $PBdim == 2 ]]
    then
	 sed "s/\(export ionicstr*=\).*\$/\1$ionicstr/
             s/\(export epsin*=\).*\$/\1$epsin/" \
		 "$DelphiDir"/DELPHI_Memb.pbp > DELPHI.pbp
    fi

}

make_sites () 
{

    "$CpHDIR"/scripts/make_sites $1 TMP_protein.pdb > ./${runname}.sites
    #
    # Correct TMP_CpHMD.top with NTR and CTR
    #
    if egrep NT ${runname}.sites; then 

        awk '
        BEGIN{
            top=ARGV[1]; 
            while (getline < top)
            if ($5=="H3") ntr[++j]=$3;
            close(top)}; 
        END{
            while (getline < top) {
            write=1; 
              for (i in ntr) {
                if ($3==ntr[i] && $5=="N") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
                if ($3==ntr[i] && $5=="H1") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
                if ($3==ntr[i] && $5=="H2") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
                if ($3==ntr[i] && $5=="H3") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
                if ($3==ntr[i] && $5=="CA") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
                if ($3==ntr[i] && $5=="C") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
                if ($3==ntr[i] && $5=="O") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
              };
              if (write) print $0}
        }' TMP_CpHMD.top > TMP_aux.top
    else
        mv TMP_CpHMD.top TMP_aux.top
    fi
    if egrep CT ${runname}.sites; then 
        awk '
        BEGIN{
            top=ARGV[1]; 
            while (getline < top)
            if ($5=="O1") ctr[++j]=$3;
            close(top)}; 
        END{
            while (getline < top) {
            write=1; 
              for (i in ctr) {
                if ($3==ctr[i] && $5=="C") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="O1") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="O2") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HO11") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HO12") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HO21") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HO22") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
              };
              if (write) print $0}
        }' TMP_aux.top > TMP_CpHMD_charge.top
    else
        mv TMP_aux.top TMP_CpHMD_charge.top
    fi
    
}

make_delphi_DB()
{
    echo "atom__res_radius_" > ./DataBaseT.siz
    echo "atom__resnumbc_charge_" > ./DataBaseT.crg
    for type in crg siz ; do
	## Insert the terminals in the database
	awk  '$2~"CT" {print}' $DelphiDir/DataBaseT_toolarge.${type} >> ./DataBaseT.${type}
	awk  '$2~"NT" {print}' $DelphiDir/DataBaseT_toolarge.${type} >> ./DataBaseT.${type}
	##
	## Insert other residues
	for db_res in `awk '/ATOM/ {print substr($0,18,2)}' TMP_solute.pdb | sort | uniq`;do	
	    if [[ ! -z `awk -v d=$db_res '$2~d {print}' "$DelphiDir"/DataBaseT_toolarge.${type}` ]] ; then
		awk -v d=$db_res '$2~d {print}' $DelphiDir/DataBaseT_toolarge.${type} >> ./DataBaseT.${type}
	    else
		message E  "Residue identifier $db_res is not on Delphi Database_toolarge. Program will crash"
	    fi	
	done
    done

    ## make CRG_FILE ##
    awk 'NF==3 {printf"%-6s%-9s%6.3f\n", $1,$2,0.000}' ./DataBaseT.crg > CRG_FILE
    
    # Generate charges database for Delphi
    "$DelphiDir"/gen_charge.awk ${runname}.sites \
        TMP_CpHMD_charge.top TMP_delphi.crg
}


check_nchains () #SCAMPOS: new function to account for several chains #15-11-2011
{
    
    nchains=`awk 'BEGIN{chain_old=" "};/ATOM/{chain_now=substr($0,22,1); if (chain_now != chain_old) {chain_old=chain_now; if (chain_now != " ") n++}}END{print n}' TMP_protein.pdb`
    if [[ $nchains -eq 0 ]]; then
        message W  "No chain indicator found in file $PDBin. One single chain considered."
    else
        message W  "Number of chains found by reading chain indicator in file $PDBin is $nchains."    
    fi
}

set_termini()
{

    #################################################################
    #SCAMPOS: New ter reading to account for several chains
    #14-11-2011
    if [[ $nchains -gt 1 ]]; then
        #Assess position of Nter and Cter on sites list 
        NTpos=( `awk '/NT/{print NR}' ${runname}.sites` )
        CTpos=( `awk '/CT/{print NR}' ${runname}.sites` )

        #Assess position of titrating Nter and titrating Cter 
        #on Nter list and Cter list
	# log removed a print last_c on "; chain_old=chain_now}}END"
	
        NTres=`awk 'BEGIN{chain_old=" "};/ATOM/{chain_now=substr($0,22,1); if (chain_now != chain_old) {new_chain=1; chain_old=chain_now}; if (chain_now != " ") if (new_chain && $3=="N"){new_chain=0; print substr($0,23,4)}}' TMP_protein.pdb`
	CTres=`awk 'BEGIN{chain_old=" "};/ATOM/{if ($3 == "C") last_c=substr($0,23,4); chain_now=substr($0,22,1); if (chain_now != chain_old) {if (chain_old != " "); chain_old=chain_now}}END{if (chain_old != " ") print last_c}' TMP_protein.pdb`	      
       
        NTindex=( `awk -v res="$NTres" 'BEGIN{split(res,r," ")};/NT/{for (i in r) if ($1==r[i]) {print i-1; break}}' ${runname}.sites` )
        CTindex=( `awk -v res="$CTres" 'BEGIN{split(res,r," ")};/CT/{for (i in r) if ($1==r[i]) {print i-1; break}}' ${runname}.sites` )

        message W  "Number of Nterminus titrating is: ${#NTindex[@]}."
        message W  "Number of Cterminus titrating is: ${#CTindex[@]}."

        #Check if Nterminus and Cterminus are given correctly in .pHmdp SC 30-11-2011
        nNtit=$(($nchains-${#NTindex[@]}))
        nCtit=$(($nchains-${#CTindex[@]}))

        message W  "Number of Nterminus not titrating is: $nNtit."
        message W  "Number of Cterminus not titrating is: $nCtit."

        if [[ $nNtit -ne 0 ]]; then
        #Attribute state to the Nter and Cter not titrating
            n=0
            for ((i=0; i<$nchains; i++))
            do
                Titrating="FALSE"
                if [[ ${NTindex[@]} != "" ]]; then
                    for j in ${NTindex[@]}
                    do
                        if [[ $i -eq $j ]]; then
                            Titrating="TRUE"
                            break
                        fi
                    done
                fi
         # This is FF hack specific. It is now working with modified GMX 4.0.7
         # Capping Conditions
                if [[ $Titrating == "FALSE" ]]; then
                #No .pHmdp: export Nterminus=(CAP CHG CHG) (e. g.)
                    case ${Nterminus[$n]} in
                        CAP )           NTer[$i]=6 ;;
                        CAPpro )        NTer[$i]=11 ;; 
                        CHG )           NTer[$i]=3 ;; 
                        NEU | TAU1 )    NTer[$i]=0 ;;
                        TAU2 )          NTer[$i]=1 ;;
                        TAU3 | TAU4 )    NTer[$i]=2 ;;
                        * ) message E "Nterminus = ${Nterminus[$n]} is not valid. Check .pHmdp file" ;;
                    esac
                    ((++n))
                fi
            done
        fi

        if [[ $nCtit -ne 0 ]]; then

#            if [[ $nCtit -ne ${#Cterminus[@]} ]]; then
#                message E  "Cterminus array contains wrong number of elements (${#Cterminus[@]}). By inspecting file ${runname}.sites, $nCtit number of elements was expected. Please check ${blockname}.pHmdp."
#            fi

            n=0
            for ((i=0; i<$nchains; i++))
            do 
                Titrating="FALSE"
                if [[ ${CTindex[@]} != "" ]]; then
                    for j in ${CTindex[@]}
                    do
                        if [[ $i -eq $j ]]; then
                            Titrating="TRUE"
                            break
                        fi
                    done
                fi
         # This is FF hack specific. It is now working with modified GMX 5.1.5
         # Capping Conditions
                if [[ $Titrating == "FALSE" ]]; then
                    case ${Cterminus[$n]} in 
                        CAP | CAPpro )  CTer[$i]=7 ;; 
                        CHG )           CTer[$i]=4 ;; 
                        NEU | TAU1 )    CTer[$i]=0 ;;
                        TAU2 )          CTer[$i]=1 ;;
                        TAU3 )          CTer[$i]=2 ;;
                        TAU4 )          CTer[$i]=3 ;;
                        REGC )          CTer[$i]=6 ;;
                        REGN )          CTer[$i]=5 ;;
                        * ) message E "Cterminus = ${Cterminus[$n]} is not valid. Check .pHmdp file" ;;
                    esac
                    ((++n))
                fi
            done
        fi
#################################################################
    else #If number of chains is 1
        
    # This is FF hack specific. It is now working with modified GMX 4.0.7
    # Capping Conditions
        if ! egrep NT ${runname}.sites; then 
            if [[ ${#Nterminus[@]} -ne 1 ]]; then #SC 30-11-2011
                message E  "Nterminus array contains wrong number of elements (${#Nterminus[@]}). One element was expected. Please check ${blockname}.pHmdp."   
            fi
            case $Nterminus in 
                CAP )           nter=6 ;;
                CAPpro )        nter=11 ;; 
                CHG )           nter=3 ;; 
                NEU | TAU1 )    nter=0 ;;
                TAU2 )          nter=1 ;;
                TAU3 | TAU4 )    nter=2 ;;
                * ) message E "Nterminus = $Nterminus is not valid. Check .pHmdp file" ;;
            esac
        fi
        if ! egrep CT ${runname}.sites; then 
            if [[ ${#Cterminus[@]} -ne 1 ]]; then #SC 30-11-2011
                message E  "Cterminus array contains wrong number of elements (${#Cterminus[@]}). One element was expected. Please check ${blockname}.pHmdp."
            fi  
            case $Cterminus in 
                CAP | CAPpro )  cter=7 ;; 
                CHG )           cter=4 ;; 
                NEU | TAU1 )    cter=0 ;;
                TAU2 )          cter=1 ;;
                TAU3 )          cter=2 ;;
                TAU4 )          cter=3 ;;
                * ) echo message E "Cterminus = $Cterminus is not valid. Check .pHmdp file" ;;
            esac
        fi
fi
}

run_PBMC()
{
        if [ ${PBdim} -eq 2 ]; then
        # Multi-step centering procedure:
        # 1- center with the tail of lipid 1
            echo -e "Onetail\nSolute" | "$GroDIR" trjconv \
                -f TMP_effective.gro -s TMP_CpHMD.tpr \
                -o TMP_effective_aux1.gro -n TMP_CpHMD.ndx \
                -center -pbc atom
        # 
        # 2- center with all tails of the monolayer that has lipid 1
            echo -e "Monotail\nSolute" | "$GroDIR" trjconv \
                -f TMP_effective_aux1.gro -s TMP_CpHMD.tpr \
                -o TMP_effective_aux2.gro -n TMP_CpHMD.ndx \
                -center -pbc atom
        # 
        # 3- center with all tails of all lipids in the system
            echo -e "Bitail\nSolute" | "$GroDIR" trjconv \
                -f TMP_effective_aux2.gro -s TMP_CpHMD.tpr \
                -o TMP_effective_aux3.gro -n TMP_CpHMD.ndx \
                -center -pbc atom
        #
        # 4- center with all lipids in the system
            echo -e "Solute\nSolute" | "$GroDIR" trjconv \
                -f TMP_effective_aux3.gro -s TMP_CpHMD.tpr \
                -o TMP_${runname}.gro -n TMP_CpHMD.ndx \
                -center -pbc atom
        else
	    if [[ $nchains -gt 1 || $fixgro -eq 1 ]] ; then 
		#
		# In case multiple objects, use fixgro!
		# A fixgro variable was created on pHmdp for rare cases
		# with multiple objects but a single chain
		echo -e "Solute\nSolute" | "$GroDIR" trjconv \
		     -f TMP_effective.gro -s TMP_CpHMD.tpr \
		     -o TMP_effective_aux1.gro -n TMP_CpHMD.ndx \
		     -pbc mol -center
		#
		"$CpHDIR"/scripts/fixgro TMP_effective_aux1.gro \
			 ${runname}.fixgro > TMP_${runname}.gro
	    else
		# Making protein whole removing PBC:
		echo "Protein" | "$GroDIR" trjconv -f TMP_effective.gro \
					   -o TMP_${runname}.gro -s TMP_CpHMD.tpr \
					   -n TMP_CpHMD.ndx -pbc mol -quiet
            fi
            #
	    #DEBUG - store PDB trajectory format for later analysis
        
            "$GroDIR" editconf -f TMP_${runname}.gro \
                   -o TMP_aux_debug.pdb
            cat TMP_aux_debug.pdb >> ${blockname}_PQR.pdb
        
	fi
    #
    # Run delphiT and Delphi: # MM 09/01/2012
        "$DelphiDir"/delphiT $nCPU ${runname} ${PBdim} \
		    >TMP_delphi.out 2>TMP_delphi.err

    ### Check for fortran STOP within the output ###
	if grep -q "FORTRAN STOP" TMP_delphi.err
	then
	    message E "Fortran Stop on delphi step, something went wrong! on PB/MC cycle."
	fi
    # Run cconvert
    "$DelphiDir"/convert ${runname}.pkcrg ${runname}.g $temp > ${runname}.dat
    #
    # Run PETIT
    newseed=$((seed + $(                                          
    if test -f $runname.mocc ;                                                                               
    then                                                                                             
      wc -l $runname.mocc | awk '{print $1}'                                                                 
    else                                                                                             
      echo 0                                                                                         
    fi                                                                                               
  )                                                                                                  
))  
    echo "$pH,$pH,1" -E "$pot,$pot,1" -T $temp -c 2 -r $newseed -q 1000 100000
    "$PetitDIR"/petit -H "$pH,$pH,1" -E "$pot,$pot,1" -T $temp \
        -c 2 -r $newseed -q 1000 100000 <${runname}.dat \
        >TMP_MCarlo.out 2>TMP_MCarlo.err 
    #
    # Removing PB related auxiliary files
    rm -f ${runname}.{summ,pqr*,dat,pkcrg,g,pkint,out} \
        *.potat TMP_aux* TMP_mead_out TMP_mead_error 

}

write_fractions ()
{
    # Write the occupation files
    awk '
    BEGIN{
    
      # Read petit output (CE/MC):
      while (getline < "TMP_MCarlo.out")
      {
        # Read occ (from CE/MC):
        if ($0 ~ /^f/)
        {
          nsites = NF - 1 ;
          for(i = 2 ; i <= NF ; i++) tocc[i-1] = occ[i-1] = $i ;
        }
        # Read mocc (from CE/MC):
        if ($0 ~ /^\./ && $0 !~ /tot/) tmocc[$4+1] = mocc[$4+1] = $5 ;
      }
      close("TMP_MCarlo.out") ;
    
      # Write tocc and tmocc:
      for (i = 1 ; i <= nsites ; i++)
      {
        printf ("%d ", tocc[i]) >> "TMP_CpHMD.occ" ;
        printf ("%f ", tmocc[i]) >> "TMP_CpHMD.mocc" ;
      }
      printf("\n") >> "TMP_CpHMD.occ" ;
      printf("\n") >> "TMP_CpHMD.mocc" ;
    
    }'

}

build_topology ()
{
    
 # Set termini if they are titrating
    
    if [[ $nchains -gt 1 ]]; then
#################################################################
#SCAMPOS: New ter reading to account for several chains
#14-11-2011
   
        n=0
        if [[ ${NTindex[@]} != "" ]]; then
            for i in ${NTindex[@]}
            do
                NTer[$i]=`tail -n -1 TMP_CpHMD.occ | awk -v i=${NTpos[$n]} '{print $i}'`
                ((++n))
            done
        fi

        n=0
        if [[ ${CTindex[@]} != "" ]]; then
            for i in ${CTindex[@]}
            do
                CTer[$i]=`tail -n -1 TMP_CpHMD.occ | awk -v i=${CTpos[$n]} '{print $i}'`
                ((++n))
            done
        fi
#################################################################
    else
        if  egrep NT ${runname}.sites; then
            nter=`tail -n -1 TMP_CpHMD.occ | awk '{printf "%u\n", $1}'`
        fi
        if  egrep CT ${runname}.sites; then 
            cter=`tail -n -1 TMP_CpHMD.occ | awk '{printf "%u\n", $NF}'`
        fi
    fi


    # Run PDB2GMX to generate the new topology
        # Define protonation states for the sites needed. The first flag of 
        # pdbswitch ($gstaut=1) defines the use of tautomers. Currently, CpHMD only works with tautomers.
        "$CpHDIR"/scripts/pdbswitch 1 TMP_MCarlo.out \
                 TMP_protein.pdb > TMP_aux1.pdb #SC 11-11-2011

    
    if [[ $nchains -gt 1 ]]; then
        #
	# Print to the .err file the N- and C-terminus protonation states used
        (for ((i=0; i<$nchains; i++)); do 
             echo -e "${NTer[$i]}\n${CTer[$i]}"; done)
	#
        (for ((i=0; i<$nchains; i++)); do 
         echo -e "${NTer[$i]}\n${CTer[$i]}\n$water"; done) | \
            "$GroDIR" pdb2gmx -f TMP_aux1.pdb \
		      -p TMP_CpHMD.top -o TMP_aux3.gro -i TMP_posre.itp \
		      -ignh -ter -ff $ffID  -water select -quiet -merge all
	#
	# Unify Protein name
        #sed -i 's/Protein_chain_A/Protein        /g' TMP_CpHMD.top
        sed -i 's/\(.*\)_chain_A/Protein        /g' TMP_CpHMD.top
    else
	#
        echo -e "$nter\n$cter\n$water"  | \
	    "$GroDIR" pdb2gmx -f TMP_aux1.pdb \
		      -p TMP_CpHMD.top -o TMP_aux3.gro -i TMP_posre.itp \
		      -ignh -ter -ff $ffID -water select -quiet -merge all
        sed -i 's/RNA/Protein        /g' TMP_CpHMD.top
    fi
    #  
    ######################################################################
    #
    #
    awk -v s="$TOPSTRING" '/moleculetype/, $0 ~ s' TMP_CpHMD.top \
        > TMP_CpHMD.toptit #SC 28-11-2011
    #
    cat TMP_CpHMD.topbegin TMP_CpHMD.toptit TMP_CpHMD.topend | \
        awk '{if($1~/Protein/){print substr($1,1,7), $2, $3}else{print $0}}' \
        > TMP_aux2.top  #SC 28-11-2011 MM 17-01-2012

    # Fix the topology with the assigned rules
    # Testing for RNA TS 20/06/2023 ; AMBER RNA FF works with different rules than GROMOS
    # Temporarily not using fix_topology until rules are added or an alternaive is achieved
    # The topology files are not empty since the bonds do not require specbonds
    mv TMP_aux2.top TMP_CpHMD.top
    #$CpHDIR/scripts/fix_topology TMP_aux2.top \
	#$RULEin > TMP_CpHMD.top 
    #
    # SC 25-11-2011
    #several rules need a solution; 
    #If RULEin is defined as an array, will it work ${RULEin[@]}? 
    #SC 13-12-2011

##     # Fix the topology using fix_dendrimer_top. LCSF 2011-10-03
##     if [[ $CpHModule == dendrimer ]]; then
##         mv TMP_CpHMD.top TMP_aux2.top
##         $CpHDIR/scripts/fix_dendrimer_top \
##             $CpHDIR/scripts/${RULEdendr} TMP_aux2.top > TMP_CpHMD.top
##     fi

    # Correct POSRE file in case of using Position Restraints
    if [[ -f $PosRe && $PosRe != *(" ") ]]; then 
        cp -f $PosRe TMP_posre.itp
    fi
    #
    # debugging the final topologies MM
    # cp TMP_CpHMD.top topology_debug_${Cycle}.top
    # Housekeeping
    rm -f TMP_aux* \#*

}

run_dynamics ()
{
    "$GroDIR" grompp -f TMP_$1.mdp -po TMP_$1_out.mdp \
        -c TMP_$2.gro -p TMP_CpHMD.top -pp TMP_processed.top \
        -n TMP_CpHMD.ndx -o TMP_$1.tpr -maxwarn 1000 -quiet

    if [ $? -eq 1 ]; then
        echo "Job cancelled due to process error."
        scancel ${SLURM_JOB_ID}
    fi
    if [ -f ${runname}_plumed.dat ] && [[ $1 == "effective" ]] && [ $Cycle -eq "1" ] && [[ "${plumed}" == "yes" ]];then
       sed '/RESTART #/d' ${runname}_plumed.dat > ${runname}_first.dat
       $mdrun -s TMP_$1.tpr -x TMP_$1.xtc -c TMP_$1.gro \
            -e TMP_$1.edr -g TMP_$1.log -o TMP_$1.trr -quiet \
            -plumed ${runname}_first.dat -nice 19 # SC-14-12-2011
        #
        if [ $? -eq 1 ]; then
            echo "Job cancelled due to process error."
            scancel ${SLURM_JOB_ID}
        fi

       
    elif [ -f ${runname}_plumed.dat ] && [[ $1 == "effective" ]] && [ $Cycle -eq "1" ] && [[ "${plumed}" == "grid" ]];then
        sed '/RESTART #/d' ${runname}_plumed.dat > ${runname}_first.dat
        sed 's/GRID_RFILE=.*//' ${runname}_first.dat  > ${runname}_final.dat

        $mdrun -s TMP_$1.tpr -x TMP_$1.xtc -c TMP_$1.gro \
            -e TMP_$1.edr -g TMP_$1.log -o TMP_$1.trr -quiet \
            -plumed ${runname}_final.dat -nice 19 # SC-14-12-2011
        #
        if [ $? -eq 1 ]; then
            echo "Job cancelled due to process error."
            scancel ${SLURM_JOB_ID}
        fi
        mv HILLS HILLS_$Cycle

    elif [ -f ${runname}_plumed.dat ] && [[ $1 == "effective" ]] && [ $Cycle -gt "1" ] && [[ "${plumed}" == "grid" ]] ;then
        $mdrun -s TMP_$1.tpr -x TMP_$1.xtc -c TMP_$1.gro \
            -e TMP_$1.edr -g TMP_$1.log -o TMP_$1.trr -quiet \
            -plumed ${runname}_plumed.dat -nice 19 # SC-14-12-2011
        #
        
        if [ $? -eq 1 ]; then
            echo "Job cancelled due to process error."
            scancel ${SLURM_JOB_ID}
        fi
        mv HILLS HILLS_$Cycle

    elif [ -f ${runname}_plumed.dat ] && [[ $1 == "effective" ]] && [ $Cycle -gt "1" ] && [[ "${plumed}" == "yes" ]]  || [[ "${plumed}" == "static" ]];then
        $mdrun -s TMP_$1.tpr -x TMP_$1.xtc -c TMP_$1.gro \
            -e TMP_$1.edr -g TMP_$1.log -o TMP_$1.trr -quiet \
            -plumed ${runname}_plumed.dat -nice 19 # SC-14-12-2011
        #

        if [ $? -eq 1 ]; then
            echo "Job cancelled due to process error."
            scancel ${SLURM_JOB_ID}
        fi

    else
        $mdrun -s TMP_$1.tpr -x TMP_$1.xtc -c TMP_$1.gro \
            -e TMP_$1.edr -g TMP_$1.log -o TMP_$1.trr -quiet \
            -nice 19 # SC-14-12-2011
        #
        if [ $? -eq 1 ]; then
            echo "Job cancelled due to process error."
            scancel ${SLURM_JOB_ID}
        fi
    fi


    rm -f \#*
}

run_relaxation () #SC 28-11-2011
{
    #Solvent relaxation
    run_dynamics relax effective

    #Prepare input GRO for dynamics
    # Solvent marker is hardcoded in this version
    SOL1st="SOL"
    awk -v s=$SOL1st '$1 ~ s {exit};{print $0}' TMP_effective.gro > TMP_aux.gro
    awk -v s=$SOL1st '$1 ~ s {a=1};a'  TMP_relax.gro >> TMP_aux.gro

#    mv -f TMP_relax.gro TMP_relax_DEBUG.gro
    mv -f TMP_aux.gro TMP_relax.gro

}


data_append ()
{
    if [ $Cycle -eq $InitCycle ]; then
        InitTime=`echo $sim_time-$WriteTime | bc -l`
        initxtc=""
        initedr=""
    else
        initxtc="TMP_CpHMD.xtc"
        initedr="TMP_CpHMD.edr"
    fi
    # Append .edr files
    echo -e "$InitTime\n`echo $sim_time-$WriteTime | bc -l`" | \
        "$GroDIR" eneconv  -o TMP_CpHMD.edr -f $initedr \
                  TMP_effective.edr -settime -quiet
    #
    # Append .xtc files
    echo -e "$InitTime\n`echo $sim_time-$WriteTime | bc -l`" | \
        "$GroDIR" trjcat -f  $initxtc TMP_effective.xtc  \
                  -o TMP_CpHMD.xtc -settime

    #
    # Append HILLS for GRID setup
    if [ -f HILLS_$Cycle ] && [[ "${plumed}" == "grid" ]]
    then
	    cat HILLS_$Cycle >> HILLS_curr_seg
    fi
    
    # Append colvar and HILLS files
    # It seems that PLUMED has no issue with the repeated timestamps between cycles.
    # The HILLS height seems to be converging to 0.
    # Also the metad.bias seems to be fine. 
    #if [ -f ${runname}_first.dat ] && [ $Cycle -eq $InitCycle ]
    #then
    #    added_time=`echo $sim_time-$WriteTime | bc -l`
    #    # In case formatting is important
    #    #awk -v t=$added_time '$1 ~ "#!" {print} ; BEGIN {OFS=""} !/#!/ {printf "%15s",$1+t; for(i=2;i<=NF;i++) printf "%20s",$i; printf "\n"}' $hills  > TMP_CpHMD_hills
    #    #awk -v t=$added_time '$1 ~ "#!" {print} ; BEGIN {OFS=""} !/#!/ {printf "%15s",$1+t; for(i=2;i<=NF;i++) printf "%20s",$i; printf "\n"}' $colvar_name > TMP_CpHMD_colvar
    #    #awk -v t=$added_time '$1 ~ "#!" {print} ; !/#!/ {$1+=t} {print $0}' $colvar_name > TMP_CpHMD_colvar
    #    #awk -v t=$added_time '$1 ~ "#!" {print} ; !/#!/ {$1+=t} {print $0}' $hills  > TMP_CpHMD_hills
    #
    #    #mv TMP_CpHMD_colvar $colvar_name
    #    #mv TMP_CpHMD_hills  $hills
    #
    #
    #elif [ -f ${runname}_plumed.dat ] && [ $Cycle -gt $InitCycle ]
    #then
    #    added_time=`echo $sim_time-$WriteTime | bc -l`
    #    # Time 0 of colvar needs to be removed. It is a repeated configuration, the same as the last one of the previous block.
    #    #awk -v t=$added_time '!/#!/ && $1!=0 && $1+=t' $colvar_name > TMP_CpHMD_colvar
    #    #awk -v t=$added_time '!/#!/ && $1+=t' $hills  > TMP_CpHMD_hills
    #    #awk -v t=$added_time '$1==0' $hills > TMP_CpHMD_hills       
    #    #awk -v t=$added_time '$1 ~ "#!" {print} ; BEGIN {OFS=""} !/#!/ {printf "%15s",$1+t; for(i=2;i<=NF;i++) printf "%20s",$i; printf "\n"}' $hills > TMP_CpHMD_hills
    #    #awk -v t=$added_time '$1 ~ "#!" {print} ; BEGIN {OFS=""} !/#!/ && $1!=0 {printf "%15s",$1+t; for(i=2;i<=NF;i++) printf "%20s",$i; printf "\n"}' $colvar_name > TMP_CpHMD_colvar
    #
    #    #cat TMP_CpHMD_colvar >> $colvar_name
    #    #cat TMP_CpHMD_hills  >>  $hills
    #
    #fi
    ##
    # Append and backup remaining files
    if [ -f TMP_effective.log -o -f TMP_effective0.log ]
    then
        cat TMP_effective*.log >> TMP_CpHMD.log
    fi
    cp -f TMP_effective.gro TMP_CpHMD.gro
    if [ -f pullx.xvg ]; then
	for e in x f ; do
	    awk -v t=`echo "$EffectiveSteps*$TimeStep*($Cycle-1)" | bc -l` \
		'!/^\#/ && !/^\@/ && $1!=0.0000{print $1+t,$2}' pull$e.xvg \
		>> TMP_CpHMD_pull$e.xvg
	done
    fi
    rm -f TMP_aux* \#* TMP_effective*.log pull?.xvg
}

clean_up ()
{
    rm -rf ${runname}.{summ,g,pkint,pkcrg,dat,pqr,out} state*\
        *.st aminoacids.dat FF.dat specbond.dat ff* traj.trr \
        ${runname}_cpu* TMP_{statepqr.in,${runname}.pqr,effective,relax,CpHMD,aux,mead,MCarlo,template,posre,processed,pqr}*
}
