#!/bin/bash
cd `dirname $0`
time ../../run -nb_threads 6 -from LONDON\ HEATHROW -to NEW\ YORK\ JOHN\ F\.\ KENNEDY -departure_time_min 10252012000000 -departure_time_max 11042012000000 -arrival_time_min 11052012000000 -arrival_time_max 11142012000000 -max_layover 34400 -vacation_time_min 700000 -vacation_time_max 904800 -vacation_airports SEATTLE NASHVILLE MEXICO\ CITY MIAMI ISTANBUL ATLANTA GUATEMALA SPOKANE\ INTL -flights flights.txt -alliances alliances.txt -work_hard_file work_hard.txt -play_hard_file play_hard.txt
