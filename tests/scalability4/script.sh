#!/bin/bash
./run -nb_threads 12 -from BOSTON -to MUNICH -departure_time_min 04102012000000 -departure_time_max 04142012000000 -arrival_time_min 05102012000000 -arrival_time_max 05142012000000 -max_layover 10000 -vacation_time_min 500000 -vacation_time_max 700000 -vacation_airports ATLANTA LONDON\ HEATHROW PHILADELPHIA RIO\ DE\ JANEIRO TORONTO $*
