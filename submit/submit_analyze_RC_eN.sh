#! /bin/env csh

set LocalPath = $PWD

cp ./macros/analyze_eD_DIS_eicrecon_RC_eN.C ./submit_workdir/.

cp ./macros/analyze_eD_DIS_eicrecon_RC_eN_C* ./submit_workdir/.

set Exec = ${LocalPath}/submit/run_analyze_submit_RC_eN.sh

####### Initialize condor file
echo ""  > ./submit/CondorFileAnalyze_RC_eN
echo "Universe     = vanilla" >> ./submit/CondorFileAnalyze_RC_eN
echo "Executable   = ${Exec}" >> ./submit/CondorFileAnalyze_RC_eN
echo "getenv = true" >> ./submit/CondorFileAnalyze_RC_eN

# Output Directory

set Output = ${LocalPath}/jobs/log
set Output_err = ${LocalPath}/jobs/err

set FileListPath = ${LocalPath}/fileLists/fileList_eD_10x130_DIS_eN_10M_events_new.list
set OutFilePath = ${LocalPath}/production/output_RC_eN_submit_work


set OutFile = ${Output}/Analyze_0_RC_eN.out
set ErrFile = ${Output_err}/Analyze_0_RC_eN.err

set Args = ( ${FileListPath} ${OutFilePath} )

echo "" >> ./submit/CondorFileAnalyze_RC_eN
echo "Output       = ${OutFile}" >> ./submit/CondorFileAnalyze_RC_eN
echo "Error        = ${ErrFile}" >> ./submit/CondorFileAnalyze_RC_eN
echo "Arguments    = ${Args}" >> ./submit/CondorFileAnalyze_RC_eN
echo "Queue" >> ./submit/CondorFileAnalyze_RC_eN   

condor_submit ./submit/CondorFileAnalyze_RC_eN


