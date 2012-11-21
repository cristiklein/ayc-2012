#!/bin/bash
cd `dirname $0`
time ../../run -nb_threads 4 -from BOSTON -to MUNICH -departure_time_min 04102012000000 -departure_time_max 04142012000000 -arrival_time_min 05102012000000 -arrival_time_max 05142012000000 -max_layover 24400 -vacation_time_min 500000 -vacation_time_max 700000 -flights flights.txt -alliances alliances.txt -work_hard_file work_hard.txt -play_hard_file play_hard.txt -vacation_airports ATLANTA LONDON\ HEATHROW PHILADELPHIA RIO\ DE\ JANEIRO TORONTO
