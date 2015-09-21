#!/bin/bash

# This demo script demonstrates projection filesystem operation
echo "Start projections filesystem demonstration"

# Folder with data
MOUNT_POINT=./data

# This folder will in future be placed in program directory not routinely seen by a user.
REMOUNT_POINT=./mnt

# Mount original data folder to program folder to enable transitive behavior (emulation of "mounting folder to itself" state)
sudo mount --bind $MOUNT_POINT $REMOUNT_POINT

# Activate FUSE filesystem on mnt folder
python3 iontorrent.py $MOUNT_POINT $REMOUNT_POINT 2>&1 > /dev/null &

echo "Press any key to terminate"

read input

echo "Ending demonstration"
sudo umount data
sudo umount mnt