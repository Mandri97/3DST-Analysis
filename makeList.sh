# This script is useful to list the files to be used in the analysis
# TODO: Remove the file with same name before running this script
# and the number of production
for i in $(seq 1 20); do
    FILE_PATH="/pnfs/dune/persistent/users/gyang/3DST/dump/standardGeo10/PROD$i/"
    ls $FILE_PATH > prod$i.txt

    # TODO: Change the output file name
    awk -v file_path="$FILE_PATH" '{print file_path$0}' prod$i.txt >> Geo10Full.txt
    awk -v file_path="$FILE_PATH" '{print file_path$0}' prod$i.txt > PROD$i.txt
    rm -rf prod$i.txt
done
