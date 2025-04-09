#!/bin/bash

R_SCRIPT="model.R"
LOG_FILE="model.log"
TIME_LOG="timestamp.log"
CLEAR_RDS="false"

for arg in "$@"; do
    case $arg in
        --script=*)
            R_SCRIPT="${arg#*=}"
            shift ;;
        --log_file=*)
            LOG_FILE="${arg#*=}"
            shift ;;
        --clear_rds=*)
            CLEAR_RDS="${arg#*=}"
            shift ;;
        *)
            echo "[ERROR] Unknown argument: $arg"
            exit 1 ;;
    esac
done

if [ ! -f "$R_SCRIPT" ]; then
  echo "[ERROR] File '$R_SCRIPT' not found"
  exit 1
fi

# clear cached rds files
shopt -s nullglob
rds_files=(*.rds)

if [[ "$CLEAR_RDS" == "true" ]]; then
  if (( ${#rds_files[@]} > 0 )); then
    rm -f "${rds_files[@]}"
    echo "[WARNING] --clear_rds=true: All cached RDS files deleted."
  else
    echo "[INFO] --clear_rds=true: No RDS files to delete."
  fi
else
  echo "[INFO] --clear_rds=false: RDS files left untouched."
fi

(
  START_TIME=$(date '+%Y-%m-%d %H:%M:%S')
  START_EPOCH=$(date '+%s')

  echo "[INFO] Job started at: $START_TIME" > "$TIME_LOG"

  nohup Rscript "$R_SCRIPT" > "$LOG_FILE" 2>&1
  EXIT_CODE=$?

  END_TIME=$(date '+%Y-%m-%d %H:%M:%S')
  END_EPOCH=$(date '+%s')
  RUN_TIME=$((END_EPOCH - START_EPOCH))

  {
    echo "[INFO] Job ended at: $END_TIME"
    echo "[INFO] Runtime: $RUN_TIME seconds"
    echo "[INFO] Exit code: $EXIT_CODE"
  } >> "$TIME_LOG"

) >/dev/null 2>&1 & 
JOB_PID=$!
disown "$JOB_PID"

echo "[INFO] PID: $JOB_PID"
