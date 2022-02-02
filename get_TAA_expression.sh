#!/bin/bash

$RES_FILE="RUN1/abundance_genelevel_TPMtrans.txt"
##
##Edit above this block
##-------------------------------------------------------------------
##Set1
grep -w "BIRC5" $RES_FILE
grep -w "CTAG1B" $RES_FILE
grep -w "PRAME" $RES_FILE
grep -w "SSX2" $RES_FILE
grep -w "WT1" $RES_FILE

##Set2
grep -w "MAGEA4" $RES_FILE
grep -w "BIRC5" $RES_FILE
grep -w "PRAME" $RES_FILE
grep -w "HPVE6" $RES_FILE
grep -w "HPVE7" $RES_FILE
