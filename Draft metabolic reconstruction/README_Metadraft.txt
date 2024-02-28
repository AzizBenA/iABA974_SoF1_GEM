1- Open the windows cmd.

2- Redirect to the folder using this command:
cd cbmpy-metadraft

3- Activate the environment:
conda activate metadraft3

4- Run the executable file with this command:
runwin.bat


Make it run on ubuntu:

1- Open Xlaunch and run an X SERVER.

2- Get your ip address with this command:

ip route | grep default | awk '{print $3}'

export DISPLAY=192.168.80.1:0

3- Repeat the same steps as windows

4- Use this command to run the GUI:
sh ./run.sh

