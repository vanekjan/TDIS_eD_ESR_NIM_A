#! /bin/env csh

set LocalPath = $PWD

cp ./macros/analyze_eD_DIS_eicrecon_RC_eP.C ./submit_workdir/.

cp ./macros/analyze_eD_DIS_eicrecon_RC_eP_C* ./submit_workdir/.

set Exec = ${LocalPath}/submit/run_analyze_submit_RC_eP.sh

####### Initialize condor file
echo ""  > ./submit/CondorFileAnalyze_RC_eP
echo "Universe     = vanilla" >> ./submit/CondorFileAnalyze_RC_eP
echo "Executable   = ${Exec}" >> ./submit/CondorFileAnalyze_RC_eP
echo "getenv = true" >> ./submit/CondorFileAnalyze_RC_eP

# Output Directory
set Output = ${LocalPath}/jobs/log
set Output_err = ${LocalPath}/jobs/err

set FileListPath = ${LocalPath}/fileLists/fileList_eD_10x130_DIS_eN_10M_events_new.list
#set FileListPath = ${LocalPath}/fileLists/fileList_eD_10x130_DIS_eP_test.list
set OutFilePath = ${LocalPath}/production/output_RC_eP_submit_work


set OutFile = ${Output}/Analyze_0_RC_eP.out
set ErrFile = ${Output_err}/Analyze_0_RC_eP.err

set Args = ( ${FileListPath} ${OutFilePath} )

echo "" >> ./submit/CondorFileAnalyze_RC_eP
echo "Output       = ${OutFile}" >> ./submit/CondorFileAnalyze_RC_eP
echo "Error        = ${ErrFile}" >> ./submit/CondorFileAnalyze_RC_eP
echo "Arguments    = ${Args}" >> ./submit/CondorFileAnalyze_RC_eP
echo "Queue" >> ./submit/CondorFileAnalyze_RC_eP   

condor_submit ./submit/CondorFileAnalyze_RC_eP

